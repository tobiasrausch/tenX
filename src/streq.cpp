/*
============================================================================
Strand-Seq Watson-Crick Counter
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


struct Config {
  bool outwc;
  bool outpre;
  bool hasRegionFile;
  bool hasPreFile;
  unsigned short minMapQual;
  uint32_t window;
  uint32_t segment;
  boost::filesystem::path watsonRatio;
  boost::filesystem::path preOutput;
  boost::filesystem::path regions;
  boost::filesystem::path loadPre;
  std::vector<boost::filesystem::path> files;
};

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
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(1000000), "window length")
    ("segment,s", boost::program_options::value<uint32_t>(&c.segment)->default_value(10000), "segment size")
    ("loadpre,l", boost::program_options::value<boost::filesystem::path>(&c.loadPre), "load preprocessing info")
    ;

  boost::program_options::options_description region("Region options");
  region.add_options()
    ("regions,r", boost::program_options::value<boost::filesystem::path>(&c.regions), "bed file with regions to analyze")
    ;

  boost::program_options::options_description output("Output options");
  output.add_options()
    ("preout,p", boost::program_options::value<boost::filesystem::path>(&c.preOutput), "output preprocessing info")
    ("wcratio,a", boost::program_options::value<boost::filesystem::path>(&c.watsonRatio), "output file for WC ratio")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input bam file")
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
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <strand.seq.bam>" << std::endl;
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

  // Check preprocessing file
  if (vm.count("loadpre")) {
    if (!(boost::filesystem::exists(c.loadPre) && boost::filesystem::is_regular_file(c.loadPre) && boost::filesystem::file_size(c.loadPre))) {
      std::cerr << "Pre-processing information is missing: " << c.loadPre.string() << std::endl;
      return 1;
    }
    c.hasPreFile = true;
  } else c.hasPreFile = false;

  // Output WC ratio
  if (vm.count("wcratio")) c.outwc = true;
  else c.outwc = false;
  if (vm.count("preout")) c.outpre = true;
  else c.outpre = false;

  // Load bam file
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);


  // Calculate valid windows
  typedef uint32_t TPos;
  typedef uint32_t TFileIndex;
  typedef std::pair<TPos, TFileIndex> TPosFilePair;
  typedef std::vector<TPosFilePair> TWindowList;
  typedef std::vector<TWindowList> TGenomicWindows;
  TGenomicWindows watsonWindows;
  watsonWindows.resize(hdr->n_targets);
  TGenomicWindows crickWindows;
  crickWindows.resize(hdr->n_targets);
  TGenomicWindows wcWindows;
  wcWindows.resize(hdr->n_targets);
  if (!c.hasPreFile) {
    // Parse bam (contig by contig)
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      std::cout << "Preprocessing " << c.files[file_c].string() << std::endl;
      for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
	if (hdr->target_len[refIndex] < c.window) continue;
	uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
	typedef std::vector<uint32_t> TCounter;
	TCounter watsonCount;
	TCounter crickCount;
	watsonCount.resize(bins, 0);
	crickCount.resize(bins, 0);
	
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  
	  int32_t pos = rec->core.pos + halfAlignmentLength(rec);
	  if (rec->core.flag & BAM_FREAD1) 
	    if (rec->core.flag & BAM_FREVERSE) ++crickCount[(int) (pos / c.window)];
	    else ++watsonCount[(int) (pos / c.window)];
	  else
	    if (rec->core.flag & BAM_FREVERSE) ++watsonCount[(int) (pos / c.window)];
	    else ++crickCount[(int) (pos / c.window)];
	}
	
	// Get reads per window
	TCounter support;
	support.resize(bins, 0);
	TCounter::iterator itSupport = support.begin();
	TCounter::const_iterator itWatson = watsonCount.begin();
	TCounter::const_iterator itCrick =crickCount.begin();
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++itSupport, ++bin) 
	  *itSupport = *itWatson + *itCrick;
	std::sort(support.begin(), support.end());
	uint32_t lowerCutoff = support[(int) (bins/100)];
	uint32_t upperCutoff = support[(int) (99*bins/100)];
	support.clear();
      
	// Get Watson Ratio
	itWatson = watsonCount.begin();
	itCrick =crickCount.begin();
	typedef std::vector<double> TRatio;
	TRatio wRatio;
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  // At least 1 read every 10,000 bases
	  if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) wRatio.push_back(((double) *itWatson / (double) (sup)));
	}
	std::sort(wRatio.begin(), wRatio.end());
	uint32_t wRatioSize = wRatio.size();
	if (!wRatioSize) continue;
	
	// Categorize chromosomes (exclude chromosomes with recombination events)
	double lowerW = wRatio[(int) (wRatioSize/10)];
	double upperW = wRatio[(int) (9*wRatioSize/10)];
	itWatson = watsonCount.begin();
	itCrick =crickCount.begin();
	if ((lowerW > 0.8) && (upperW > 0.8)) {
	  // Hom. Watson
	  for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	    uint32_t sup = *itWatson + *itCrick;
	    if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	      double wTmpRatio = (double) *itWatson / (double) (sup);
	      if (wTmpRatio > 0.8) watsonWindows[refIndex].push_back(std::make_pair(bin, file_c));
	    }
	  }
	} else if ((lowerW < 0.2) && (upperW < 0.2)) {
	  // Hom. Crick
	  for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	    uint32_t sup = *itWatson + *itCrick;
	    if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	      double wTmpRatio = (double) *itWatson / (double) (sup);
	      if (wTmpRatio < 0.2) crickWindows[refIndex].push_back(std::make_pair(bin, file_c));
	    }
	  }
	} else if ((lowerW > 0.3) && (upperW < 0.7)) {
	  //WC
	  for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	    uint32_t sup = *itWatson + *itCrick;
	    if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	      double wTmpRatio = (double) *itWatson / (double) (sup);
	      if ((wTmpRatio > 0.3) && (wTmpRatio < 0.7)) wcWindows[refIndex].push_back(std::make_pair(bin, file_c));
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }    
  } else {
    std::ifstream preFile(c.loadPre.string().c_str(), std::ifstream::in);
    if (preFile.is_open()) {
      while (preFile.good()) {
	std::string line;
	getline(preFile, line);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(line, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  int32_t ind = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t refIndex = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t bin = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t fileInd = boost::lexical_cast<int32_t>(*tokIter++);
	  if (ind == 0) watsonWindows[refIndex].push_back(std::make_pair(bin, fileInd));
	  else if (ind == 1) crickWindows[refIndex].push_back(std::make_pair(bin, fileInd));
	  else if (ind == 2) wcWindows[refIndex].push_back(std::make_pair(bin, fileInd));
	}
      }
      preFile.close();
    }
  }

  // Sort windows
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    std::sort(watsonWindows[refIndex].begin(), watsonWindows[refIndex].end());
    std::sort(crickWindows[refIndex].begin(), crickWindows[refIndex].end());
    std::sort(wcWindows[refIndex].begin(), wcWindows[refIndex].end());
  }

  // Dump pre-processing information
  if (c.outpre) {
    std::ofstream ofile(c.preOutput.string().c_str());
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      for(TWindowList::iterator iW = watsonWindows[refIndex].begin(); iW != watsonWindows[refIndex].end(); ++iW) ofile << "0\t" << refIndex << '\t' << iW->first << '\t' << iW->second << '\t' << std::endl;
      for(TWindowList::iterator iC = crickWindows[refIndex].begin(); iC != crickWindows[refIndex].end(); ++iC) ofile << "1\t" << refIndex << '\t' << iC->first << '\t' << iC->second << '\t' << std::endl;
      for(TWindowList::iterator iWC = wcWindows[refIndex].begin(); iWC != wcWindows[refIndex].end(); ++iWC) ofile << "2\t" << refIndex << '\t' << iWC->first << '\t' << iWC->second << '\t' << std::endl;
    }
    ofile.close();
  }

  // Store intervals
  typedef std::pair<uint32_t, uint32_t> TInterval;
  typedef std::map<TInterval, std::string> TIntervalMap;
  typedef std::vector<TIntervalMap> TGenomicIntervals;
  TGenomicIntervals genomicIntervals;
  genomicIntervals.resize(hdr->n_targets);
  if (c.hasRegionFile) {
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
	    uint32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	    uint32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	    std::string intId = *tokIter++;
	    genomicIntervals[tid].insert(std::make_pair(std::make_pair(start, end), intId));
	  }
	}
      }
      regionFile.close();
    }
  } else {
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      std::string intId(hdr->target_name[refIndex]);
      intId += ":" + boost::lexical_cast<std::string>(0) + "-" + boost::lexical_cast<std::string>(hdr->target_len[refIndex]);
      genomicIntervals[refIndex].insert(std::make_pair(std::make_pair(0, hdr->target_len[refIndex]), intId));
    }
  }

  // Output watson ratios
  if (c.outwc) {
    std::ofstream ofile(c.watsonRatio.string().c_str());
    ofile << "chr\tstart\tend\twratio\tsupport\ttype\tid" << std::endl;
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      for(TIntervalMap::const_iterator itR = genomicIntervals[refIndex].begin(); itR != genomicIntervals[refIndex].end(); ++itR) {
	std::string intervalName = itR->second;
	uint32_t intervalSize = itR->first.second - itR->first.first;
	uint32_t indS = (int) (itR->first.first / c.window);
	uint32_t bins = intervalSize / c.segment + 1;
	for(std::size_t i = 0; i<3; ++i) {
	  typedef std::vector<uint32_t> TCounter;
	  TCounter watsonCount;
	  TCounter crickCount;
	  watsonCount.resize(bins, 0);
	  crickCount.resize(bins, 0);

	  TWindowList::const_iterator itWindows;
	  TWindowList::const_iterator itWindowsEnd;
	  if (i == 0) {
	    itWindows = std::lower_bound(watsonWindows[refIndex].begin(), watsonWindows[refIndex].end(), std::make_pair(indS, (uint32_t) 0));
	    itWindowsEnd = watsonWindows[refIndex].end();
	  } else if (i == 1) {
	    itWindows = std::lower_bound(crickWindows[refIndex].begin(), crickWindows[refIndex].end(), std::make_pair(indS, (uint32_t) 0));
	    itWindowsEnd = crickWindows[refIndex].end();
	  } else if (i == 2) {
	    itWindows = std::lower_bound(wcWindows[refIndex].begin(), wcWindows[refIndex].end(), std::make_pair(indS, (uint32_t) 0));
	    itWindowsEnd = wcWindows[refIndex].end();
	  }
	  for(;itWindows != itWindowsEnd; ++itWindows) {
	    if (itWindows->first * c.window >= itR->first.second) break;
	    int32_t regionStart = std::max(itWindows->first * c.window, itR->first.first);
	    int32_t regionEnd = std::min((itWindows->first + 1) * c.window, itR->first.second);
	    hts_itr_t* iter = sam_itr_queryi(idx[itWindows->second], refIndex, regionStart, regionEnd);
	    bam1_t* rec = bam_init1();
	    while (sam_itr_next(samfile[itWindows->second], iter, rec) >= 0) {
	      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	      
	      int32_t pos = rec->core.pos + halfAlignmentLength(rec);
	      if ((pos >= (int32_t) itR->first.first) && (pos < (int32_t) itR->first.second)) {
		uint32_t binInd = (uint32_t) ((pos - itR->first.first) / c.segment);
		if (rec->core.flag & BAM_FREAD1) 
		  if (rec->core.flag & BAM_FREVERSE) ++crickCount[binInd];
		  else ++watsonCount[binInd];
		else
		  if (rec->core.flag & BAM_FREVERSE) ++watsonCount[binInd];
		  else ++crickCount[binInd];
	      }
	    }
	  }

	  TCounter::const_iterator itWatson = watsonCount.begin();
	  TCounter::const_iterator itCrick =crickCount.begin();
	  for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	    uint32_t sup = *itWatson + *itCrick;
	    if (sup > 0) {
	      double wRatio = ((double) *itWatson / (double) (sup));
	      int32_t regionStart = bin * c.segment + itR->first.first;
	      int32_t regionEnd = std::min((uint32_t) ((bin + 1) * c.segment + itR->first.first), itR->first.second);
	      ofile << hdr->target_name[refIndex] << '\t' << regionStart << '\t' << regionEnd << '\t' << wRatio << '\t' << sup << '\t';
	      if (i == 0) ofile << "Watson" << '\t';
	      else if (i == 1) ofile << "Crick" << '\t';
	      else if (i == 2) ofile << "WatsonCrick" << '\t';
	      ofile << intervalName << std::endl;
	    }
	  }
	}
      }
    }
    ofile.close();
  }

  // Close bam
  bam_hdr_destroy(hdr);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }

  return 0;
}
