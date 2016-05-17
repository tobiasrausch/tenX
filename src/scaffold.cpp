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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/vcf.h>
#include <htslib/sam.h>

#include <stdio.h>


#define BARCODE_LENGTH 14

struct Config {
  bool outbar;
  bool outccomp;
  bool outspan;
  bool hasRegionFile;
  unsigned short minMapQual;
  uint32_t contiglen;
  uint32_t window;
  uint32_t linkedreads;
  uint32_t minsupport;
  uint32_t stat;
  boost::filesystem::path bamfile;
  boost::filesystem::path regions;
  boost::filesystem::path barcodeTable;
  boost::filesystem::path componentTable;
  boost::filesystem::path spanningTree;
};

struct VertexProp {
  int32_t chr;
};

struct SInterval {
  uint32_t chr;
  uint32_t start;
  uint32_t end;

  SInterval(uint32_t c, uint32_t s, uint32_t e) : chr(c), start(s), end(e) {}
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
  for(;itC != itEnd; ++itC, ++pos) hash += (uint32_t) (*itC * std::pow((double) 4, (int32_t) pos));
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


template<typename TConfig, typename TBarcodeChrCount>
inline void _countContigHits(TConfig const& c, samFile* samfile, hts_idx_t* idx, bam_hdr_t* hdr, TBarcodeChrCount& barChrCount) {
  typedef typename TBarcodeChrCount::key_type TBarcodeChr;

  // Parse bam (contig by contig)
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

	uint8_t *miptr = bam_aux_get(rec, "MI");
	if (miptr) {
	  uint8_t *bxptr = bam_aux_get(rec, "BX");
	  if (bxptr) {
	    char* bx = (char*) (bxptr + 1);
	    // Get barcode id
	    uint32_t hash = _encodeBarcode(boost::to_upper_copy(std::string(bx)));

	    // Get chromosome id
	    int32_t chrId = rec->core.tid + 1;
	    if (rec->core.pos + halfAlignmentLength(rec) < midpoint) chrId = -1 * chrId;

	    // Insert the barcode-chromosome count
	    TBarcodeChr bc = std::make_pair(hash, chrId);
	    typename TBarcodeChrCount::iterator itBC = barChrCount.find(bc);
	    if (itBC != barChrCount.end()) itBC->second += 1;
	    else barChrCount[bc] = 1;
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
  }
}


template<typename TConfig, typename TBarcodeChrCount, typename TIntervalList>
inline void _countRegionHits(TConfig const& c, samFile* samfile, hts_idx_t* idx, bam_hdr_t* hdr, TBarcodeChrCount& barChrCount, TIntervalList& genomeIntervals) {
  typedef typename TBarcodeChrCount::key_type TBarcodeChr;
  typedef typename TIntervalList::value_type TSInterval;

  // Parse intervals
  typedef boost::icl::interval_set<uint32_t> TChrIntervals;
  typedef typename TChrIntervals::interval_type TIVal;
  typedef std::vector<TChrIntervals> TRegionsGenome;
  TRegionsGenome regions;
  regions.resize(hdr->n_targets);
  std::ifstream regionFile(c.regions.string().c_str(), std::ifstream::in);
  if (regionFile.is_open()) {
    while (regionFile.good()) {
      std::string line;
      getline(regionFile, line);
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t,;");
      Tokenizer tokens(line, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	std::string chrName = *tokIter++;
	int32_t tid = bam_name2id(hdr, chrName.c_str());
	if (tid >= 0) {
	  int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	  if ((end - start) < (int32_t) c.contiglen) continue;
	  regions[tid].insert(TIVal::right_open(start, end));
	}
      }
    }
    regionFile.close();
  }


  // Parse bam (contig by contig)
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    for(TChrIntervals::const_iterator itR = regions[refIndex].begin(); itR != regions[refIndex].end(); ++itR) {
      uint32_t intervalSize = itR->upper() - itR->lower();
      uint32_t intervalindex = genomeIntervals.size();
      genomeIntervals.push_back(TSInterval(refIndex, itR->lower(), itR->upper()));
      int32_t midpoint = (int32_t) (intervalSize / 2) + itR->lower();
      for (int contigLoc = 0; contigLoc < 2; ++contigLoc) {
	int32_t regionStart = itR->lower();
	int32_t regionEnd = std::min(midpoint, (int32_t) (itR->lower() + c.window));
	if (contigLoc) {
	  regionStart = std::max(midpoint + 1, (int32_t) itR->upper() - (int32_t) c.window);
	  regionEnd = std::min(itR->upper(), hdr->target_len[refIndex]);
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
	    int32_t chrId = intervalindex + 1;
	    if (rec->core.pos + halfAlignmentLength(rec) < midpoint) chrId = -1 * chrId;

	    // Insert the barcode-chromosome count
	    TBarcodeChr bc = std::make_pair(hash, chrId);
	    typename TBarcodeChrCount::iterator itBC = barChrCount.find(bc);
	    if (itBC != barChrCount.end()) itBC->second += 1;
	    else barChrCount[bc] = 1;
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }
  }
}


int main(int argc, char **argv) {
  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("contiglen,c", boost::program_options::value<uint32_t>(&c.contiglen)->default_value(4000), "min. contig length")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(25000), "max. contig start/end search window (molecule length)")
    ("linkedreads,l", boost::program_options::value<uint32_t>(&c.linkedreads)->default_value(3), "min. linked reads per contig")
    ("minsupport,m", boost::program_options::value<uint32_t>(&c.minsupport)->default_value(15), "min. #barcode contig-contig support")
    ;
  boost::program_options::options_description region("Region options");
  region.add_options()
    ("regions,r", boost::program_options::value<boost::filesystem::path>(&c.regions), "bed file with regions to analyze")
    ;

  boost::program_options::options_description output("Output options");
  output.add_options()
    ("barcodetable,b", boost::program_options::value<boost::filesystem::path>(&c.barcodeTable), "output barcode,#contigs,#reads count table")
    ("concomponents,p", boost::program_options::value<boost::filesystem::path>(&c.componentTable), "output connected components")
    ("spanningtree,s", boost::program_options::value<boost::filesystem::path>(&c.spanningTree), "output min. spanning trees")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input bam file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(region).add(output).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(region).add(output);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample.10X.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Check regions file
  if (vm.count("regions")) {
    if (!(boost::filesystem::exists(c.regions) && boost::filesystem::is_regular_file(c.regions) && boost::filesystem::file_size(c.regions))) {
      std::cerr << "Region file is missing: " << c.regions.string() << std::endl;
      return 1;
    }
    c.hasRegionFile = true;
  } else c.hasRegionFile = false;

  // Output barcode table
  if (vm.count("barcodetable")) c.outbar = true;
  else c.outbar = false;
  if (vm.count("concomponents")) c.outccomp = true;
  else c.outccomp = false;
  if (vm.count("spanningtree")) c.outspan = true;
  else c.outspan = false;

  // Load bam file
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  // Count barcodes per contig/region
  typedef uint32_t TBarcode;
  typedef int32_t TChr;
  typedef std::pair<TBarcode, TChr> TBarcodeChr;
  typedef boost::unordered_map<TBarcodeChr, uint32_t> TBarcodeChrCount;
  TBarcodeChrCount barChrCount;
  typedef std::vector<SInterval> TIntervalList;
  TIntervalList genomeIntervals;
  if (c.hasRegionFile) _countRegionHits(c, samfile, idx, hdr, barChrCount, genomeIntervals);
  else _countContigHits(c, samfile, idx, hdr, barChrCount);

  // Count barcodes
  typedef boost::unordered_map<TBarcode, std::size_t> TBarcodeIndex;
  TBarcodeIndex barIndex;
  std::size_t barCount = 0;
  for(TBarcodeChrCount::const_iterator itBC = barChrCount.begin();itBC != barChrCount.end(); ++itBC)
    if (barIndex.find(itBC->first.first) == barIndex.end()) barIndex.insert(std::make_pair(itBC->first.first, barCount++));

  // Estimate barcode,#contigs,#reads distributions
  std::vector<std::size_t> barNumContigs;
  barNumContigs.resize(barCount, 0);
  std::vector<std::size_t> barNumReads;
  barNumReads.resize(barCount, 0);
  for(TBarcodeChrCount::const_iterator itBC = barChrCount.begin();itBC != barChrCount.end(); ++itBC) {
    TBarcodeIndex::iterator itBarIndex = barIndex.find(itBC->first.first);
    if ((itBC->first.second >= 0) || (barChrCount.find(std::make_pair(itBC->first.first, std::abs(itBC->first.second))) == barChrCount.end())) ++barNumContigs[itBarIndex->second];
    barNumReads[itBarIndex->second] += itBC->second;
  }

  // Output barcode table
  if (c.outbar) {
    std::ofstream ofile(c.barcodeTable.string().c_str());
    for(TBarcodeIndex::const_iterator itBarIndex = barIndex.begin(); itBarIndex != barIndex.end(); ++itBarIndex) {
      ofile << _decodeBarcode(itBarIndex->first) << "\t" << barNumContigs[itBarIndex->second] << "\t" << barNumReads[itBarIndex->second] << std::endl;
    }
    ofile.close();
  }

  // Identify valid, spanning barcodes
  boost::unordered_set<TBarcode> barValid;
  for(TBarcodeIndex::const_iterator itBarIndex = barIndex.begin(); itBarIndex != barIndex.end(); ++itBarIndex) {
    if ((barNumContigs[itBarIndex->second] > 1) && (barNumReads[itBarIndex->second]/barNumContigs[itBarIndex->second] >= c.linkedreads)) {
      barValid.insert(itBarIndex->first);
    }
  }
  barIndex.clear();
  barNumContigs.clear();
  barNumReads.clear();

  // Summarize chromosome counts per valid barcode
  typedef std::pair<TChr, uint32_t> TChrCount;
  typedef std::vector<TChrCount> TChrCountVec;
  typedef boost::unordered_map<TBarcode, TChrCountVec> TBarcodeGenome;
  TBarcodeGenome barGenome;
  for(TBarcodeChrCount::const_iterator itBC = barChrCount.begin();itBC != barChrCount.end(); ++itBC) {
    // Valid barcode ?
    if (barValid.find(itBC->first.first) != barValid.end()) {
      bool validContig = false;
      if (itBC->second >= c.linkedreads) validContig = true;
      else {
	TBarcodeChrCount::iterator itOtherPart = barChrCount.find(std::make_pair(itBC->first.first, -1 * itBC->first.second));
	if ((itOtherPart != barChrCount.end()) && (itBC->second + itOtherPart->second >= c.linkedreads)) validContig = true;
      }
      if (validContig) {
	TBarcodeGenome::iterator itBarGenome = barGenome.find(itBC->first.first);
	if (itBarGenome == barGenome.end()) itBarGenome = barGenome.insert(std::make_pair(itBC->first.first, TChrCountVec())).first;
	itBarGenome->second.push_back(std::make_pair(itBC->first.second, itBC->second));
      }
    }
  }
  barChrCount.clear();

  // Pre-filter chromosomal links
  typedef boost::unordered_map<std::size_t, uint8_t> TChrPairFilter;
  TChrPairFilter chrPairFilter;
  for(TBarcodeGenome::const_iterator itBG = barGenome.begin(); itBG != barGenome.end(); ++itBG) {
    for(TChrCountVec::const_iterator itCVIt = itBG->second.begin(); itCVIt != itBG->second.end(); ++itCVIt) {
      TChrCountVec::const_iterator itCVItNext = itCVIt;
      ++itCVItNext;
      for(;itCVItNext != itBG->second.end(); ++itCVItNext) {
	TChr chr1 = std::min(std::abs(itCVIt->first), std::abs(itCVItNext->first));
	TChr chr2 = std::max(std::abs(itCVIt->first), std::abs(itCVItNext->first));
	if (chr1 != chr2) {
	  std::size_t pairHash = hashChrPair(std::make_pair(chr1, chr2));
	  if (chrPairFilter.find(pairHash) == chrPairFilter.end()) chrPairFilter[pairHash] = 0;
	  if (chrPairFilter[pairHash]<255) chrPairFilter[pairHash] += 1;
	}
      }
    }
  }

  // Summarize chromosomal links
  typedef std::pair<TChr, TChr> TChrPair;
  typedef boost::unordered_map<TChrPair, uint32_t> TChrPairLinks;
  TChrPairLinks chrLinks;
  for(TBarcodeGenome::const_iterator itBG = barGenome.begin();itBG != barGenome.end(); ++itBG) {
    for(TChrCountVec::const_iterator itCVIt = itBG->second.begin(); itCVIt != itBG->second.end(); ++itCVIt) {
      TChrCountVec::const_iterator itCVItNext = itCVIt;
      ++itCVItNext;
      for(;itCVItNext != itBG->second.end(); ++itCVItNext) {
	TChr chr1 = std::min(std::abs(itCVIt->first), std::abs(itCVItNext->first));
	TChr chr2 = std::max(std::abs(itCVIt->first), std::abs(itCVItNext->first));
	if (chr1 != chr2) {
	  TChrPair cP = std::make_pair(chr1, chr2);
	  std::size_t pairHash = hashChrPair(cP);
	  if (chrPairFilter[pairHash] >= c.minsupport) {
	    TChrPairLinks::iterator itChrLinks = chrLinks.find(cP);
	    if (itChrLinks == chrLinks.end()) itChrLinks = chrLinks.insert(std::make_pair(cP, 0)).first;
	    itChrLinks->second += 1;
	  }
	}
      }
    }
  }
  chrPairFilter.clear();

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
    }
  }
  chrLinks.clear();

  // Get connected components IDs
  typedef boost::unordered_set<uint32_t> TComponentSet;
  TComponentSet compIds;
  typedef std::pair<uint32_t, TChr> TCompChrPair;
  typedef std::vector<TCompChrPair> TCompChr;
  TCompChr compChr;
  for(std::size_t i = 0; i < comp.size(); ++i) {
    if (comp[i]>0) {
      compIds.insert(comp[i]);
      compChr.push_back(std::make_pair(comp[i], i));
    }
  }

  // Sort by component
  std::sort(compChr.begin(), compChr.end());

  // Output connected components
  if (c.outccomp) {
    std::ofstream ofile(c.componentTable.string().c_str());
    uint32_t lastCompID = 0;
    TCompChr::const_iterator itConnChr = compChr.begin();
    TCompChr::const_iterator itConnChrEnd = compChr.end();
    for(;itConnChr != itConnChrEnd; ++itConnChr) {
      std::string chrName;
      if (c.hasRegionFile) {
	chrName += hdr->target_name[genomeIntervals[itConnChr->second - 1].chr];
	chrName += ":";
	chrName += boost::lexical_cast<std::string>(genomeIntervals[itConnChr->second - 1].start) + "-";
	chrName += boost::lexical_cast<std::string>(genomeIntervals[itConnChr->second - 1].end);
      } else chrName += hdr->target_name[itConnChr->second - 1];
      if (itConnChr->first != lastCompID) {
	if (lastCompID) ofile << std::endl;
	ofile << "scaffold" << itConnChr->first << ": " << chrName;
	lastCompID = itConnChr->first;
      } else {
	ofile << "," << chrName;
      }
    }
    ofile << std::endl;
    ofile.close();
  }

  // Iterate connected components
  std::ofstream fout;
  if (c.outspan) fout.open(c.spanningTree.string().c_str());
  for(TComponentSet::const_iterator ici = compIds.begin(); ici != compIds.end(); ++ici) {
    typedef int32_t TEdgeWeight;
    TEdgeWeight defaultInnerContigEdge = -10000;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProp, boost::property<boost::edge_weight_t, TEdgeWeight> > Graph;
    Graph g;

    // Edge property map
    typedef boost::property_map<Graph, boost::edge_weight_t>::type TEdgeMap;
    TEdgeMap weightMap = get(boost::edge_weight, g);

    // Define the reverse map
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::unordered_map<TChr, Vertex> TChrVertexMap;
    TChrVertexMap chrVertex;

    // Create all vertices
    for(TCompChr::const_iterator itConnChr = compChr.begin();itConnChr != compChr.end(); ++itConnChr) {
      if (itConnChr->first == *ici) {
	TChrVertexMap::iterator pos;
	bool inserted;

	// Add Left Node
	Vertex u;
	boost::tie(pos, inserted) = chrVertex.insert(std::make_pair(-1 * itConnChr->second, Vertex()));
	if (inserted) {
	  u = add_vertex(g);
	  pos->second = u;
	  g[u].chr = -1 * itConnChr->second;
	} else {
	  u = pos->second;
	}

	// Add Right Node
	Vertex v;
	boost::tie(pos, inserted) = chrVertex.insert(std::make_pair(itConnChr->second, Vertex()));
	if (inserted) {
	  v = add_vertex(g);
	  pos->second = v;
	  g[v].chr = itConnChr->second;
	} else {
	  v = pos->second;
	}

	// Add inner contig edge
	boost::graph_traits<Graph>::edge_descriptor e;
	bool edgeins;
	tie(e, edgeins) = add_edge(u, v, g);
	if (edgeins) weightMap[e] = defaultInnerContigEdge;
      }
    }

    // Create edges
    for(TBarcodeGenome::const_iterator itBG = barGenome.begin();itBG != barGenome.end(); ++itBG) {
      for(TChrCountVec::const_iterator itCVIt = itBG->second.begin(); itCVIt != itBG->second.end(); ++itCVIt) {
	TChrCountVec::const_iterator itCVItNext = itCVIt;
	++itCVItNext;
	for(;itCVItNext != itBG->second.end(); ++itCVItNext) {
	  TChr chr1 = std::min(itCVIt->first, itCVItNext->first);
	  TChr chr2 = std::max(itCVIt->first, itCVItNext->first);
	  if ((comp[std::abs(chr1)] == *ici) and (comp[std::abs(chr1)] == comp[std::abs(chr2)])) {
	    if (std::abs(chr1) != std::abs(chr2)) {
	      TChrVertexMap::iterator pos;
	      bool inserted;
	      boost::tie(pos, inserted) = chrVertex.insert(std::make_pair(chr1, Vertex()));
	      Vertex v = pos->second;
	      boost::tie(pos, inserted) = chrVertex.insert(std::make_pair(chr2, Vertex()));
	      Vertex u = pos->second;
	      boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
	      boost::tie(ei, ei_end) = out_edges(v, g);
	      bool foundEdge = false;
	      for(; ei!=ei_end; ++ei) {
		if (target(*ei, g) == u) {
		  foundEdge = true;
		  //weightMap[*ei] -= std::min(itCVIt->second, itCVItNext->second);
		  weightMap[*ei] -= 1;
		  break;
		}
	      }
	      if (!foundEdge) {
		boost::graph_traits<Graph>::edge_descriptor e;
		tie(e, inserted) = add_edge(v, u, g);
		//if (inserted) weightMap[e] = -1 * std::min(itCVIt->second, itCVItNext->second);
		if (inserted) weightMap[e] = -1;
	      }
	    }
	  }
	}
      }
    }

    // Compute minimum spanning tree
    std::vector<boost::graph_traits<Graph>::edge_descriptor> treeEdges;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(treeEdges));

    // Write Graph
    if (c.outspan) {
      //write_graphviz(fout, g, boost::make_label_writer(boost::get(&VertexProp::name, g)), boost::make_label_writer(boost::get(boost::edge_weight, g)));
      // ... or ...
      fout << "graph scaffold" << *ici << " {\n" << " node[shape=\"box\"]\n";
      boost::graph_traits<Graph>::vertex_iterator viter, viter_end;
      for (boost::tie(viter, viter_end) = vertices(g); viter != viter_end; ++viter) {
	std::string chrName;
	if (c.hasRegionFile) {
	  chrName += hdr->target_name[genomeIntervals[std::abs(g[*viter].chr) - 1].chr];
	  chrName += ":";
	  chrName += boost::lexical_cast<std::string>(genomeIntervals[std::abs(g[*viter].chr) - 1].start) + "-";
	  chrName += boost::lexical_cast<std::string>(genomeIntervals[std::abs(g[*viter].chr) - 1].end);
	} else chrName += hdr->target_name[std::abs(g[*viter].chr) - 1];
	if (g[*viter].chr < 0) chrName.append("/L");
	else chrName.append("/R");
	if (*viter % 2 == 0) {
	  if (*viter == 0) fout << "subgraph _" << *viter << " {\nrank=\"min\";\n";
	  else if (*viter + 2 == num_vertices(g)) fout << "subgraph _" << *viter << " {\nrank=\"max\";\n";
	  else fout << "subgraph _" << *viter << " {\nrank=\"same\";\n";
	}
	fout << *viter << "[color=\"black\", label=\"" << chrName << "\"];\n";
	if (*viter % 2 == 1) fout << "}\n";
      }
      boost::graph_traits<Graph>::edge_iterator eiter, eiter_end;
      for (boost::tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter) {
	fout << source(*eiter, g) << " -- " << target(*eiter, g);
	if (std::find(treeEdges.begin(), treeEdges.end(), *eiter) != treeEdges.end()) fout << "[color=\"black\"";
	else fout << "[color=\"gray\"";
	fout << ", label=\"" << -1 * get(boost::edge_weight, g, *eiter) << "\", weight=" << -1 * get(boost::edge_weight, g, *eiter) << ", penwidth=" << 1.0 + ((double) (-1 * get(boost::edge_weight, g, *eiter)) / 10000.0) << "];\n";
      }
      fout << " }\n\n";
    }
  }
  if (c.outspan) fout.close();

  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  return 0;
}
