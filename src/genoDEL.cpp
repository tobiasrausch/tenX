/*
============================================================================
SV FDR
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

#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h> 
#include <stdio.h>

inline bool splitPoint(bam1_t* rec, int32_t& clipSize, int32_t& split, unsigned int qualCut) {
  // Check for single soft-clip
  unsigned int numSoftClip = 0;
  uint32_t* cigar = bam_get_cigar(rec);
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i)
    if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) ++numSoftClip;
  if (numSoftClip != 1) return false;

  // Get quality vector
  typedef std::vector<uint8_t> TQuality;
  TQuality quality;
  quality.resize(rec->core.l_qseq);
  uint8_t* qualptr = bam_get_qual(rec);
  for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
  
  // Get soft-clips
  uint32_t alen = 0;
  uint32_t lastIns = 0;
  unsigned int meanQuality = 0;
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
      alen += bam_cigar_oplen(cigar[i]) + lastIns;
      lastIns = 0;
    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
      lastIns = bam_cigar_oplen(cigar[i]);   // Only add if followed by 'M'
    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
      clipSize = bam_cigar_oplen(cigar[i]);
      split = rec->core.pos + alen;
      unsigned int qualSum = 0;
      for(unsigned int i = alen; i < (alen+clipSize); ++i) qualSum += quality[i];
      meanQuality = qualSum / clipSize;
    }
  }
  return (meanQuality>=qualCut);
}


inline void addBpCounts(bam1_t* rec, int32_t regionStart, int32_t regionEnd, std::vector<uint32_t>& bp) {
  if (rec->core.pos >= regionEnd) return;
  int32_t bpPos = rec->core.pos;
  uint32_t* cigar = bam_get_cigar(rec);
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
    int op = bam_cigar_op(cigar[i]);
    int ol = bam_cigar_oplen(cigar[i]);
    if (op == BAM_CMATCH)
      for(int k = 0; k<ol; ++k, ++bpPos) {
	if ((bpPos>=regionStart) && (bpPos<regionEnd)) ++bp[bpPos - regionStart];
      }
    else if ((op == BAM_CREF_SKIP) || (op == BAM_CDEL)) bpPos += ol;
  }
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " NA12878 <in.vcf.gz> <sample.bam>" << std::endl;
    return 1; 
  }

  // Parameter
  std::string sampleName = std::string(argv[1]);
  int32_t readLen = 150;
  int32_t bpwindow = 25;
  double covThreshold = 2.0;

  // Load bcf file
  htsFile* ifile = bcf_open(argv[2], "r");
  if (!ifile) {
    std::cerr << "Fail to load " << argv[2] << "!" << std::endl;
    return 1;
  }
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  bcf1_t* rec = bcf_init();

  // Load bam file
  samFile* samfile = sam_open(argv[3], "r");
  hts_idx_t* idx = sam_index_load(samfile, argv[3]);
  bam_hdr_t* hd = sam_hdr_read(samfile);

  // Parse VCF
  int32_t nsvend = 0;
  int32_t* svend = NULL;
  int32_t nsvt = 0;
  char* svt = NULL;
  int ngt = 0;
  int32_t* gt = NULL;
  int sampleIndex = 0;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
    if (hdr->samples[i] == sampleName) {
      sampleIndex = i;
      break;
    }
  }
  std::cout << "chr\tstart\tend\tid\tsize\tsvtype\thaplotype\tgenotype\tphasedblockid\tsrH1\tcovH1\tsrH2\tcovH2\tcalledhaplotype\tcalledgenotype" << std::endl;
  while (bcf_read(ifile, hdr, rec) == 0) {
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
	bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	if (std::string(svt) == "DEL") {
	  bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);

	  // Phased block identifier
	  std::string uniquePS = "None";
	  
	  // SV region
	  std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	  if (chrName.substr(0,3) != "chr") chrName = std::string("chr").append(chrName);
	  if ((chrName == "chrX") || (chrName == "chrY")) continue;   // Only autosomes
	    
	  int32_t regionChr = bam_name2id(hd, chrName.c_str());
	  int32_t regionStart = std::max(0, rec->pos - readLen);
	  int32_t regionEnd = (*svend) + readLen;
	  int32_t regionSize = regionEnd - regionStart;

	  // Split-read information
	  typedef std::vector<uint32_t> TSplitPoints;
	  TSplitPoints sppH1;
	  sppH1.resize(regionSize, 0);
	  TSplitPoints sppH2;
	  sppH2.resize(regionSize, 0);

	  // Base pair counts
	  typedef std::vector<uint32_t> TBPCounts;
	  TBPCounts bpH1;
	  bpH1.resize(regionSize, 0);
	  TBPCounts bpH2;
	  bpH2.resize(regionSize, 0);
	  std::size_t bpCountH1 = 0;
	  std::size_t bpCountH2 = 0;

	  // Initial breakpoints
	  int32_t svStart = rec->pos - regionStart;
	  int32_t svEnd = (*svend) - regionStart;
	  
	  hts_itr_t* iter = sam_itr_queryi(idx, regionChr, regionStart, regionEnd);
	  bam1_t* r = bam_init1();
	  bool singlePhasedBlock = true;
	  bool psSet = false;
	  while (sam_itr_next(samfile, iter, r) >= 0) {
	    if (r->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if (r->core.pos < regionStart) continue;
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

		// Store base-count
		if (hap==1) addBpCounts(r, regionStart, regionEnd, bpH1);
		else addBpCounts(r, regionStart, regionEnd, bpH2);
		
		// Store split-read
		int32_t split = -1;
		int32_t clipSize = -1;
		if (splitPoint(r, clipSize, split, 20)) {
		  if ((split >= regionStart) && (split < regionEnd)) {
		    split -= regionStart;
		    if (clipSize > (int) (log10(r->core.l_qseq) * 10)) {
		      if (hap==1) ++sppH1[split];
		      else ++sppH2[split];
		    }
		  }
		}
	      }
	    }
	  }
	  bam_destroy1(r);
	  hts_itr_destroy(iter);

	  if (singlePhasedBlock) {
	    // Search refined breakpoint
	    int32_t svstartbeg = std::max(0, (rec->pos - regionStart) - bpwindow);
	    int32_t svstartend = (rec->pos - regionStart) + bpwindow;
	    int32_t svendbeg = ((*svend) - regionStart) - bpwindow;
	    int32_t svendend = ((*svend) - regionStart) + bpwindow;
	    if (svstartend > svendbeg) {
	      int32_t mid = (svendbeg + svstartend) / 2;
	      svstartend = mid - 1;
	      svendbeg = mid;
	    }

	    // Which haplotype to search for split reads?
	    for(std::size_t i = svStart + 1; i < (std::size_t) svEnd; ++i) {
	      bpCountH1 += bpH1[i];
	      bpCountH2 += bpH2[i];
	    }
	    TSplitPoints spp = sppH1;
	    if (bpCountH1 > bpCountH2) spp = sppH2;
	    uint32_t sC = spp[svStart];
	    for(std::size_t i = svstartbeg; i < (std::size_t) svstartend; ++i) {
	      if (spp[i] > sC) {
		svStart = i;
		sC = spp[i];
	      }
	    }
	    sC = spp[svEnd];
	    for(std::size_t i = svendbeg; i < (std::size_t) svendend; ++i) {
	      if (spp[i] > sC) {
		svEnd = i;
		sC = spp[i];
	      }
	    }

	    // Re-calculate bp counts based on refined SV breakpoints
	    bpCountH1 = 0;
	    bpCountH2 = 0;
	    for(std::size_t i = svStart + 1; i < (std::size_t) svEnd; ++i) {
	      bpCountH1 += bpH1[i];
	      bpCountH2 += bpH2[i];
	    }
	  } else uniquePS = "None";
	  // Calculate coverage
	  std::size_t svSize = svEnd - svStart;
	  double covH1 = (double) bpCountH1 / (double) svSize;
	  double covH2 = (double) bpCountH2 / (double) svSize;

	  // SR support
	  uint32_t srH1 = sppH1[svStart] + sppH1[svEnd];
	  uint32_t srH2 = sppH2[svStart] + sppH2[svEnd];

	  // Call genotype
	  std::string gtstr = "None";
	  int gtcalled = -1;
	  if (uniquePS != "None") {
	    if (covH1 < covThreshold) {
	      if (covH2 < covThreshold) {
		if ((srH1 >= 2) && (srH2 >= 2)) {
		  gtstr = "1|1";
		  gtcalled = 2;
		}
	      } else {
		if ((srH1 >= 2) && (srH2 == 0)) { 
		    gtstr = "1|0";
		    gtcalled = 1;
		}
	      }
	    } else {
	      if (covH2 < covThreshold) {
		if ((srH1 == 0) && (srH2 >= 2)) {
		  gtstr = "0|1";
		  gtcalled = 1;
		}
	      } else {
		if ((srH1 == 0) && (srH2 == 0)) {
		  gtstr = "0|0";
		  gtcalled = 0;
		}
	      }
	    }
	  }
	  std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << *svend << "\t" << rec->d.id << "\t" << ((*svend) - rec->pos) << "\t" << svt << "\t" << gtval << "\t" << gt_type << "\t" << uniquePS << "\t" << srH1 << "\t" << covH1 << "\t" << srH2 << "\t" << covH2  << "\t" << gtstr << "\t" << gtcalled << std::endl;
	}
      }
    }
  }
  // Close bam
  hts_idx_destroy(idx);
  sam_close(samfile);
  bam_hdr_destroy(hd);

  // Clean-up
  free(svend);
  free(svt);
  free(gt);

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}
