/*
============================================================================
Bam SV aligner
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

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <tags.h>
#include <util.h>
#include <msa.h>
#include <split.h>
#include <junction.h>
#include <math.h> 
#include <stdio.h>

struct Config {
  int32_t minimumFlankSize;
  uint32_t mincov;
  uint32_t mapqual;
  double scoring;
  std::string sample;
  boost::filesystem::path vcffile;
  boost::filesystem::path genome;
  boost::filesystem::path bamfile;
};


int main(int argc, char **argv) {
  Config c;

  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "Sample")
    ("flanking,f", boost::program_options::value<int32_t>(&c.minimumFlankSize)->default_value(13), "breakpoint padding")
    ("mincov,m", boost::program_options::value<uint32_t>(&c.mincov)->default_value(5), "min. haploid coverage")
    ("qual,q", boost::program_options::value<uint32_t>(&c.mapqual)->default_value(20), "min. mapping quality")
    ("scoring,r", boost::program_options::value<double>(&c.scoring)->default_value(0.9), "min. carrier alignment score")
    ("vcf,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile)->default_value("sample.bcf"), "input bcf file")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <sample.bam> ..." << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Load bcf file
  htsFile* ifile = hts_open(c.vcffile.string().c_str(), "r");
  hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Load bam file
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  bam_hdr_t* hd = sam_hdr_read(samfile);

  // Parse VCF
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t ninslen = 0;
  int32_t* inslen = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int32_t ncons = 0;
  char* cons = NULL;
  int ngt = 0;
  int32_t* gt = NULL;
  int sampleIndex = 0;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
    if (hdr->samples[i] == c.sample) {
      sampleIndex = i;
      break;
    }
  }

  // Header
  std::cout << "chr\tstart\tend\tid\tinslen\thaplotype\tgenotype\tphasedblockid\tscoreH1\tscoreH2\tcalledhaplotype\tcalledgenotype" << std::endl;

  // Parse genome
  kseq_t *seq;
  int l;
  gzFile fp = gzopen(c.genome.string().c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seqname(seq->name.s);
    int32_t chrid = bcf_hdr_name2id(hdr, seqname.c_str());
    if (chrid == -1) {
      if (seqname.substr(0,3) != "chr") seqname = std::string("chr").append(seqname);
      if ((seqname == "chrX") || (seqname == "chrY")) continue;   // Only autosomes
      chrid = bcf_hdr_name2id(hdr, seqname.c_str());
      if (chrid == -1) continue;
    }
    hts_itr_t* itervcf = bcf_itr_queryi(bcfidx, chrid, 0, seq->seq.l);
    bcf1_t* rec = bcf_init1();
    while (bcf_itr_next(ifile, itervcf, rec) >= 0) {
      bcf_unpack(rec, BCF_UN_ALL);
      bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
      if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1)) {
	int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
	std::string gtval;
	if (bcf_gt_is_phased(gt[sampleIndex*2 + 1])) {
	  std::ostringstream s;
	  s << bcf_gt_allele(gt[sampleIndex*2]) << '|' << bcf_gt_allele(gt[sampleIndex*2 + 1]);
	  gtval = s.str();
	} else {
	  std::ostringstream s;
	  s << bcf_gt_allele(gt[sampleIndex*2]) << '/' << bcf_gt_allele(gt[sampleIndex*2 + 1]);
	  gtval = s.str();
	}
	if (gt_type != -1) {
	  bool precise = false;
	  if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
	  if (precise) {
	    bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	    std::string svtype(svt);
	    bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend)
;	    bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen);
	    bcf_get_info_string(hdr, rec, "CONSENSUS", &cons, &ncons);
	    std::string consensus = boost::to_upper_copy(std::string(cons));

	    // Get reference junk
	    std::string ref = boost::to_upper_copy(std::string(seq->seq.s + std::max((int32_t) 0, (int32_t) (rec->pos - consensus.size())), seq->seq.s + (*svend) + consensus.size()));

	    // Find breakpoint
	    typedef boost::multi_array<char, 2> TAlign;
	    typedef typename TAlign::index TAIndex;
	    TAlign alignRef;
	    torali::AlignConfig<true, false> semiglobal;
	    torali::DnaScore<int> sc(5, -4, -5 * c.minimumFlankSize, 0);
	    torali::gotoh(consensus, ref, alignRef, semiglobal, sc);
	    TAIndex cStart, cEnd, rStart, rEnd;
	    torali::_findSplit(alignRef, cStart, cEnd, rStart, rEnd);

	    // Get only the inserted sequence
	    consensus = consensus.substr(cStart, (cEnd - cStart));

	    // Re-align the inserted sequence to the local reference to ensure it's unique
	    TAlign localAlign;
	    double locsc = torali::gotoh(consensus, ref, localAlign, semiglobal);
	    locsc /= (consensus.size() * 5);

	    // Unique inserted sequence?
	    if (locsc >= 0.5 ) continue;

	    // SV region
	    std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	    if (chrName.substr(0,3) != "chr") chrName = std::string("chr").append(chrName);
	    if ((chrName == "chrX") || (chrName == "chrY")) continue;   // Only autosomes
	    
	    // Alignment scores
	    typedef std::vector<double> TScore;
	    TScore scoreH1;
	    TScore scoreH2;
	      
	    // Scan BAM
	    int32_t regionChr = bam_name2id(hd, chrName.c_str());
	    int32_t regionStart = std::max(0, rec->pos - c.minimumFlankSize);
	    int32_t regionEnd = (*svend) + c.minimumFlankSize;
	    hts_itr_t* iter = sam_itr_queryi(idx, regionChr, regionStart, regionEnd);
	    bam1_t* r = bam_init1();

	    std::string uniquePS = "None";
	    bool singlePhasedBlock = true;
	    bool psSet = false;
	    while (sam_itr_next(samfile, iter, r) >= 0) {
	      if (r->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	      if (((r->core.pos > regionStart) || ((r->core.pos + r->core.l_qseq) < (rec->pos + c.minimumFlankSize))) && ((r->core.pos > ((*svend) - c.minimumFlankSize)) || ((r->core.pos + r->core.l_qseq) < regionEnd))) continue;
	      if (r->core.qual < c.mapqual) continue;
	      uint8_t *psptr = bam_aux_get(r, "PS");
	      if (psptr) {
		// Only consider reads belonging to a phased block
		char* ps = (char*) (psptr + 1);
		std::string rPS = std::string(ps);

		// Make sure the entire SV is spanned by one PS
		if (psSet) {
		  if (rPS != uniquePS) {
		    singlePhasedBlock = false;
		    break;
		  }
		} else {
		  psSet = true;
		  uniquePS = rPS;
		}

		// Get the haplotype
		uint8_t* hpptr = bam_aux_get(r, "HP");
		if (hpptr) {
		  int hap = bam_aux2i(hpptr);
	      
		  // Get sequence
		  std::string sequence;
		  sequence.resize(r->core.l_qseq);
		  uint8_t* seqptr = bam_get_seq(r);
		  for (int i = 0; i < r->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	      
		  // Get quality
		  uint32_t qualSum = 0;
		  uint8_t* qualptr = bam_get_qual(r);
		  for (int i = 0; i < r->core.l_qseq; ++i) qualSum += qualptr[i];
		  if ((qualSum / r->core.l_qseq) < c.mapqual) continue;
	      
		  // Compute alignment to consensus
		  TAlign align;
		  torali::AlignConfig<true, true> gapfree;
		  double sc = torali::gotoh(consensus, sequence, align, gapfree);
		  sc /= (consensus.size() * 5);
		  if (hap == 1) scoreH1.push_back(sc);
		  else scoreH2.push_back(sc);
		}
	      }
	    }
	    bam_destroy1(r);
	    hts_itr_destroy(iter);

	    if (!singlePhasedBlock) uniquePS = "None";

	    // Evaluate alignment scores
	    std::string gtstr = "None";
	    int32_t gtcalled = -1;
	    double medH1 = -1;
	    double medH2 = -1;
	    if ((scoreH1.size() >= c.mincov) && (scoreH2.size() >= c.mincov)) {
	      std::sort(scoreH1.begin(), scoreH1.end());
	      std::sort(scoreH2.begin(), scoreH2.end());
	      medH1 = scoreH1[scoreH1.size()/2];
	      medH2 = scoreH2[scoreH2.size()/2];
	      if (medH1 < c.scoring) {
		if (medH2 < c.scoring) {
		  gtstr = "0|0";
		  gtcalled = 0;
		} else {
		  gtstr = "0|1";
		  gtcalled = 1;
		}
	      } else {
		if (medH2 < c.scoring) {
		  gtstr = "1|0";
		  gtcalled = 1;
		} else {
		  gtstr = "1|1";
		  gtcalled = 2;
		}
	      }
	    }

	    // Output genotype
	    std::cout << bcf_hdr_id2name(hdr, rec->rid) << '\t' << (rec->pos + 1) << '\t' << (*svend) << '\t' << rec->d.id << '\t' << (*inslen) << '\t' << gtval << '\t' << gt_type << '\t' << uniquePS << "\t" << medH1 << "\t" << medH2 << "\t" << gtstr << "\t" << gtcalled << std::endl;
	  }	  
	}
      }
    }
    bcf_destroy(rec);
    hts_itr_destroy(itervcf);
  }
  kseq_destroy(seq);
  gzclose(fp);

  // Close bam
  hts_idx_destroy(idx);
  sam_close(samfile);
  bam_hdr_destroy(hd);

  // Clean-up
  free(svend);
  free(inslen);
  free(svt);
  free(gt);

  // Close VCF
  bcf_hdr_destroy(hdr);
  hts_idx_destroy(bcfidx);
  bcf_close(ifile);

  return 0;
}
