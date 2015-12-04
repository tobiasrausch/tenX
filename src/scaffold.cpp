/*
============================================================================
10X Deletion Genotyper
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
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h> 
#include <stdio.h>

struct Config {
  unsigned short minMapQual;
  uint32_t contiglen;
  uint32_t moleculelen;
  uint32_t linkedreads;
  uint32_t stat;
  boost::filesystem::path bamfile;
};

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
    ("contiglen,c", boost::program_options::value<uint32_t>(&c.contiglen)->default_value(500), "min. contig length")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("moleculelen,m", boost::program_options::value<uint32_t>(&c.moleculelen)->default_value(75000), "max. molecule length")
    ("linkedreads,l", boost::program_options::value<uint32_t>(&c.linkedreads)->default_value(3), "min. linked reads per chr")
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
      int32_t regionEnd = std::min(midpoint, (int32_t) c.moleculelen);
      if (contigLoc) {
	regionStart = std::max(midpoint, (int32_t) hdr->target_len[refIndex] - (int32_t) c.moleculelen);
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
	  TSpanningBarcode::iterator sBIt = spanBar.find(hash);
	  if (sBIt == spanBar.end()) {
	    spanBar.insert(std::make_pair(hash, false)).first;
	    lastBarChr.insert(std::make_pair(hash, rec->core.tid));
	  } else if (!sBIt->second) {
	    if (rec->core.tid != lastBarChr.find(hash)->second) sBIt->second = true;
	  }
	  
	  // Insert the barcode-chromosome count
	  TBarcodeChr bc = std::make_pair(hash, chrId);
	  TBarcodeChrCount::iterator itBC = barChrCount.find(bc);
	  if (itBC != barChrCount.end()) itBC->second += 1;
	  else barChrCount.insert(std::make_pair(bc, 1));
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
  TBarcodeChrCount::const_iterator itBC = barChrCount.begin();
  TBarcodeChrCount::const_iterator itBCEnd = barChrCount.end();
  for(;itBC != itBCEnd; ++itBC) {
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

  // Summarize chromosomal links
  typedef std::pair<TChr, TChr> TChrPair;
  typedef boost::unordered_map<TChrPair, uint32_t> TChrPairLinks;
  TChrPairLinks chrLinks;
  TBarcodeGenome::const_iterator itBG = barGenome.begin();
  TBarcodeGenome::const_iterator itBGEnd = barGenome.end();
  for(;itBG != itBGEnd; ++itBG) {
    TChrCountVec::const_iterator itCVIt = itBG->second.begin();
    TChrCountVec::const_iterator itCVItEnd = itBG->second.end();
    for(;itCVIt != itCVItEnd; ++itCVIt) {
      TChrCountVec::const_iterator itCVItNext = itCVIt;
      ++itCVItNext;
      for(;itCVItNext != itCVItEnd; ++itCVItNext) {
	TChrPair cP = std::make_pair(std::min(itCVIt->first, itCVItNext->first), std::max(itCVIt->first, itCVItNext->first));
	TChrPairLinks::iterator itChrLinks = chrLinks.find(cP);
	if (itChrLinks == chrLinks.end()) itChrLinks = chrLinks.insert(std::make_pair(cP, 0)).first;
	if (c.stat == 0) itChrLinks->second += 1;
	else if (c.stat == 1) itChrLinks->second += itCVIt->second;
      }
    }
  }
  barGenome.clear();
  
  // Output connections
  TChrPairLinks::iterator chrIt = chrLinks.begin();
  TChrPairLinks::iterator chrItEnd = chrLinks.end();
  for(;chrIt != chrItEnd; ++chrIt) {
    std::string chr1End;
    if (chrIt->first.first < 0) chr1End = "/L";
    else chr1End = "/R";
    std::string chr2End;
    if (chrIt->first.second < 0) chr2End = "/L";
    else chr2End = "/R";
    std::string chr1Name = hdr->target_name[std::abs(chrIt->first.first) - 1];
    chr1Name = chr1Name.append(chr1End);
    std::string chr2Name = hdr->target_name[std::abs(chrIt->first.second) - 1];
    chr2Name = chr2Name.append(chr2End);
    std::cout << chr1Name << ',' << chr2Name << '\t' << chrIt->second << std::endl;
  }
  
  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  return 0;
}
