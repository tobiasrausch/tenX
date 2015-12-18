/*
============================================================================
10X Scaffolder
============================================================================
Copyright (C) 2015 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h> 
#include <stdio.h>


#define BARCODE_LENGTH 14

struct Config {
  bool outbar;
  unsigned short minMapQual;
  uint32_t contiglen;
  uint32_t window;
  uint32_t linkedreads;
  uint32_t minsupport;
  uint32_t stat;
  boost::filesystem::path bamfile;
  boost::filesystem::path barcodeTable;
};

template<typename TPair>
inline std::size_t hashChrPair(TPair const& p) {
  std::size_t seed = 0;
  boost::hash_combine(seed, p.first);
  boost::hash_combine(seed, p.second);
  return seed;
}


inline uint32_t 
_encodeBarcode(std::string const& bar) {
  typedef std::vector<unsigned int> TCount;
  TCount nucl;
  std::string::const_iterator itA = bar.begin();
  std::string::const_iterator itE = bar.end();
  for(;itA != itE; ++itA) {
    if (*itA == 'A') nucl.push_back(0);
    else if (*itA =='C') nucl.push_back(1);
    else if (*itA =='G') nucl.push_back(2);
    else if (*itA =='T') nucl.push_back(3);
  }
  uint32_t pos = 0;
  uint32_t hash = 0;
  TCount::const_iterator itC = nucl.begin();
  TCount::const_iterator itEnd = nucl.end();
  for(;itC != itEnd; ++itC, ++pos) hash += (uint32_t) (*itC * std::pow(4, pos));
  return hash;
}

inline std::string
_decodeBarcode(uint32_t hash) {
  std::string seq;
  for (std::size_t i = 0; i < BARCODE_LENGTH; ++i) {
    if (hash % 4 == 0) seq += 'A';
    else if (hash % 4 == 1) seq += 'C';
    else if (hash % 4 == 2) seq += 'G';
    else if (hash % 4 == 3) seq += 'T';
    hash /= 4;
  }
  return seq;
}

template<typename TIterator, typename TValue>
inline void
_getMedian(TIterator begin, TIterator end, TValue& median) {
  std::nth_element(begin, begin + (end - begin) / 2, end);
  median = *(begin + (end - begin) / 2);
}

template<typename TIterator, typename TValue>
inline void
_getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) {
  std::vector<TValue> absDev;
  for(;begin<end;++begin)
    absDev.push_back(std::abs((TValue)*begin - median));
  _getMedian(absDev.begin(), absDev.end(), mad);
}


inline int32_t halfAlignmentLength(bam1_t* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  unsigned int alen = 0;
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i)
    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
  return (alen / 2);
}

int main(int argc, char **argv) {
  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("contiglen,c", boost::program_options::value<uint32_t>(&c.contiglen)->default_value(1000), "min. contig length")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(25000), "max. contig start/end search window (molecule length)")
    ("linkedreads,l", boost::program_options::value<uint32_t>(&c.linkedreads)->default_value(3), "min. linked reads per chr")
    ("minsupport,m", boost::program_options::value<uint32_t>(&c.minsupport)->default_value(5), "min. scaffold support")
    ("barcodetable,b", boost::program_options::value<boost::filesystem::path>(&c.barcodeTable), "output barcode,#contigs,#reads count table")
    ("stat,s", boost::program_options::value<uint32_t>(&c.stat)->default_value(0), "0: binary link statistic, 1: sum of shared barcodes")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input bam file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample.10X.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Output barcode table
  if (vm.count("barcodetable")) c.outbar = true;
  else c.outbar = false;

  // Load bam file
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  // Parse bam (chr by chr)
  typedef uint32_t TBarcode;
  typedef int32_t TChr;
  typedef boost::unordered_map<TBarcode, bool> TSpanningBarcode;
  TSpanningBarcode spanBar;
  typedef boost::unordered_map<TBarcode, TChr> TLastBarcodeChr;
  TLastBarcodeChr lastBarChr;
  typedef std::pair<TBarcode, TChr> TBarcodeChr;
  typedef boost::unordered_map<TBarcodeChr, uint32_t> TBarcodeChrCount;
  TBarcodeChrCount barChrCount;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    if (hdr->target_len[refIndex] < c.contiglen) continue;
    int32_t midpoint = (int32_t) hdr->target_len[refIndex] / 2;
    for (int contigLoc = 0; contigLoc < 2; ++contigLoc) {
      int32_t regionStart = 0;
      int32_t regionEnd = std::min(midpoint, (int32_t) c.window);
      if (contigLoc) {
	regionStart = std::max(midpoint, (int32_t) hdr->target_len[refIndex] - (int32_t) c.window);
	regionEnd = hdr->target_len[refIndex];
      }
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, regionStart, regionEnd);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
      
	uint8_t *bxptr = bam_aux_get(rec, "BX");
	if (bxptr) {
	  char* bx = (char*) (bxptr + 1);
	  // Get barcode id
	  uint32_t hash = _encodeBarcode(boost::to_upper_copy(std::string(bx)));
	  
	  // Get chromosome id
	  int32_t chrId = rec->core.tid + 1;
	  if (rec->core.pos + halfAlignmentLength(rec) < midpoint) chrId = -1 * chrId;
	  
	  // Is this a spanning barcode?
	  if (spanBar.find(hash) == spanBar.end()) {
	    spanBar[hash] = false;
	    lastBarChr[hash] = rec->core.tid;
	  } else {
	    if ((!spanBar[hash]) && (rec->core.tid != lastBarChr[hash])) spanBar[hash] = true;
	  }
	  
	  // Insert the barcode-chromosome count
	  TBarcodeChr bc = std::make_pair(hash, chrId);
	  TBarcodeChrCount::iterator itBC = barChrCount.find(bc);
	  if (itBC != barChrCount.end()) itBC->second += 1;
	  else barChrCount[bc] = 1;
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
  }
  lastBarChr.clear();

  // Summarize chromosome counts per barcode (that do span multiple chromosomes)
  typedef std::pair<TChr, uint32_t> TChrCount;
  typedef std::vector<TChrCount> TChrCountVec;
  typedef boost::unordered_map<TBarcode, TChrCountVec> TBarcodeGenome;
  TBarcodeGenome barGenome;
  for(TBarcodeChrCount::const_iterator itBC = barChrCount.begin();itBC != barChrCount.end(); ++itBC) {
    // Spanning barcode ?
    if (spanBar.find(itBC->first.first)->second) {
      if (itBC->second >= c.linkedreads) {
	TBarcodeGenome::iterator itBarGenome = barGenome.find(itBC->first.first);
	if (itBarGenome == barGenome.end()) itBarGenome = barGenome.insert(std::make_pair(itBC->first.first, TChrCountVec())).first;
	itBarGenome->second.push_back(std::make_pair(itBC->first.second, itBC->second));
      }
    }
  }
  spanBar.clear();
  barChrCount.clear();

  // Estimate barcode,#contigs,#reads distributions and remove highly abundant barcodes and floating oligos (<50 reads)
  std::vector<std::size_t> barNumContigs;
  std::vector<std::size_t> barNumReads;
  for(TBarcodeGenome::const_iterator itBG = barGenome.begin(); itBG != barGenome.end(); ++itBG) {
    barNumContigs.push_back(itBG->second.size());
    uint32_t readsPerBarcode = 0;
    for(TChrCountVec::const_iterator itChrCount = itBG->second.begin(); itChrCount!=itBG->second.end(); ++itChrCount) readsPerBarcode += itChrCount->second;
    barNumReads.push_back(readsPerBarcode);
  }
  int32_t contigMedian = 0;
  int32_t contigMad = 0;
  _getMedian(barNumContigs.begin(), barNumContigs.end(), contigMedian);
  _getMAD(barNumContigs.begin(), barNumContigs.end(), contigMedian, contigMad);
  uint32_t barContigCutoff = (uint32_t) (contigMedian + 5 * contigMad);
  int32_t readMedian = 0;
  int32_t readMad = 0;
  _getMedian(barNumReads.begin(), barNumReads.end(), readMedian);
  _getMAD(barNumReads.begin(), barNumReads.end(), readMedian, readMad);
  uint32_t barReadCutoff = (uint32_t) (readMedian + 5 * readMad);
  barNumContigs.clear();
  barNumReads.clear();

  // Output barcode table
  if (c.outbar) {
    std::ofstream ofile(c.barcodeTable.string().c_str());
    for(TBarcodeGenome::const_iterator itBG = barGenome.begin(); itBG != barGenome.end(); ++itBG) {
      uint32_t readsPerBarcode = 0;
      for(TChrCountVec::const_iterator itChrCount = itBG->second.begin(); itChrCount!=itBG->second.end(); ++itChrCount) readsPerBarcode += itChrCount->second;
      ofile << _decodeBarcode(itBG->first) << "\t" << itBG->second.size() << "\t" << readsPerBarcode << std::endl;
    }
    ofile.close();
  }

  // Chromosomal links pre-filter
  typedef boost::unordered_map<std::size_t, uint8_t> TChrPairFilter;
  TChrPairFilter chrPairFilter;
  typedef boost::unordered_map<TBarcode, bool> TValidBarcode;
  TValidBarcode validBarcode;
  for(TBarcodeGenome::const_iterator itBG = barGenome.begin(); itBG != barGenome.end(); ++itBG) {
    uint32_t readsPerBarcode = 0;
    for(TChrCountVec::const_iterator itChrCount = itBG->second.begin(); itChrCount!=itBG->second.end(); ++itChrCount) readsPerBarcode += itChrCount->second;
    if ((itBG->second.size() <= barContigCutoff) && (readsPerBarcode<=barReadCutoff) && (readsPerBarcode>50)) validBarcode[itBG->first] = true;
    else validBarcode[itBG->first] = false;
    if (validBarcode[itBG->first]) {
      for(TChrCountVec::const_iterator itCVIt = itBG->second.begin(); itCVIt != itBG->second.end(); ++itCVIt) {
	TChrCountVec::const_iterator itCVItNext = itCVIt;
	++itCVItNext;
	for(;itCVItNext != itBG->second.end(); ++itCVItNext) {
	  std::size_t pairHash = hashChrPair(std::make_pair(std::min(itCVIt->first, itCVItNext->first), std::max(itCVIt->first, itCVItNext->first)));
	  if (chrPairFilter.find(pairHash) == chrPairFilter.end()) chrPairFilter[pairHash] = 0;
	  if (c.stat == 0) {
	    if (chrPairFilter[pairHash]<255) chrPairFilter[pairHash] += 1;
	  } else if (c.stat == 1) {
	    uint32_t incr = std::min(itCVIt->second, itCVItNext->second);
	    if ((chrPairFilter[pairHash] + incr) < 255) chrPairFilter[pairHash] += incr;
	  }
	}
      }
    }
  }

  // Summarize chromosomal links
  typedef std::pair<TChr, TChr> TChrPair;
  typedef boost::unordered_map<TChrPair, uint32_t> TChrPairLinks;
  TChrPairLinks chrLinks;
  for(TBarcodeGenome::const_iterator itBG = barGenome.begin();itBG != barGenome.end(); ++itBG) {
    if (validBarcode[itBG->first]) {
      for(TChrCountVec::const_iterator itCVIt = itBG->second.begin(); itCVIt != itBG->second.end(); ++itCVIt) {
	TChrCountVec::const_iterator itCVItNext = itCVIt;
	++itCVItNext;
	for(;itCVItNext != itBG->second.end(); ++itCVItNext) {
	  TChrPair cP = std::make_pair(std::min(itCVIt->first, itCVItNext->first), std::max(itCVIt->first, itCVItNext->first));
	  std::size_t pairHash = hashChrPair(cP);
	  if (chrPairFilter[pairHash] >= c.minsupport) {
	    TChrPairLinks::iterator itChrLinks = chrLinks.find(cP);
	    if (itChrLinks == chrLinks.end()) itChrLinks = chrLinks.insert(std::make_pair(cP, 0)).first;
	    if (c.stat == 0) itChrLinks->second += 1;
	    else if (c.stat == 1) itChrLinks->second += std::min(itCVIt->second, itCVItNext->second);
	  }
	}
      }
    }
  }
  barGenome.clear();
  chrPairFilter.clear();
  validBarcode.clear();

  // Compute connected components
  typedef std::vector<uint32_t> TComponent;
  TComponent comp;
  comp.resize(hdr->n_targets + 1, 0);
  uint32_t numComp = 0;
  TChrPairLinks::iterator chrIt = chrLinks.begin();
  TChrPairLinks::iterator chrItEnd = chrLinks.end();
  for(;chrIt != chrItEnd; ++chrIt) {
    if (chrIt->second >= c.minsupport) {
      TChr chr1 = std::abs(chrIt->first.first);
      TChr chr2 = std::abs(chrIt->first.second);
      if (chr1 != chr2) {
	uint32_t compIndex = 0;
	if (!comp[chr1]) {
	  if (!comp[chr2]) {
	    // Both vertices have no component yet
	    compIndex = ++numComp;
	    comp[chr1] = compIndex;
	    comp[chr2] = compIndex;
	  } else {
	    compIndex = comp[chr2];
	    comp[chr1] = compIndex;
	  }
	} else {
	  if (!comp[chr2]) {
	    compIndex = comp[chr1];
	    comp[chr2] = compIndex;
	  } else {
	    // Both vertices already have a component ID
	    if (comp[chr1] == comp[chr2]) {
	      compIndex = comp[chr1];
	    } else {
	      // Merge components
	      compIndex = comp[chr1];
	      uint32_t otherIndex = comp[chr2];
	      if (otherIndex < compIndex) {
		compIndex = comp[chr2];
		otherIndex = comp[chr1];
	      }
	      // Re-label other index
	      for(std::size_t i = 0; i <= comp.size(); ++i) {
		if (otherIndex == comp[i]) comp[i] = compIndex;
	      }
	    }
	  }
	}
      }
      //std::string chr1End;
      //if (chrIt->first.first < 0) chr1End = "/L";
      //else chr1End = "/R";
      //std::string chr2End;
      //if (chrIt->first.second < 0) chr2End = "/L";
      //else chr2End = "/R";
      //std::string chr1Name = hdr->target_name[std::abs(chrIt->first.first) - 1];
      //chr1Name = chr1Name.append(chr1End);
      //std::string chr2Name = hdr->target_name[std::abs(chrIt->first.second) - 1];
      //chr2Name = chr2Name.append(chr2End);
      //std::cout << chr1Name << ',' << chr2Name << '\t' << chrIt->second << std::endl;
    }
  }

  // Summarize connected components
  typedef std::pair<uint32_t, TChr> TCompChrPair;
  typedef std::vector<TCompChrPair> TCompChr;
  TCompChr compChr;
  for(std::size_t i = 0; i < comp.size(); ++i)
    if (comp[i]>0) compChr.push_back(std::make_pair(comp[i], i));
  comp.clear();

  // Sort by component
  std::sort(compChr.begin(), compChr.end());

  // Output connected components:
  uint32_t lastCompID = 0;
  TCompChr::const_iterator itConnChr = compChr.begin();
  TCompChr::const_iterator itConnChrEnd = compChr.end();
  for(;itConnChr != itConnChrEnd; ++itConnChr) {
    if (itConnChr->first != lastCompID) {
      if (lastCompID) std::cout << std::endl;
      std::cout << "Scaffold" << itConnChr->first << ": " << hdr->target_name[itConnChr->second - 1];
      lastCompID = itConnChr->first;
    } else {
      std::cout << "," << hdr->target_name[itConnChr->second - 1];
    }
  }
  std::cout << std::endl;
  
  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  return 0;
}
