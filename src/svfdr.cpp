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
#include <set>
#include <vector>
#include <algorithm>
#include <htslib/vcf.h>
#include <htslib/sam.h>

inline bool splitPoint(bam1_t* rec, int32_t& split) {
  uint32_t alen = 0;
  uint32_t lastIns = 0;
  uint32_t numSoftClip = 0;
  uint32_t* cigar = bam_get_cigar(rec);
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i) {
    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
      alen += bam_cigar_oplen(cigar[i]) + lastIns;
      lastIns = 0;
    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
      lastIns = bam_cigar_oplen(cigar[i]);   // Only add if followed by 'M'
    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
      ++numSoftClip;
      split = rec->core.pos + alen;
    }
  }
  return (numSoftClip == 1);
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
  int32_t readLen = 100;

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
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1)) {
      int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
      if (gt_type != 0) {
	bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
	if (std::string(svt) == "DEL") {
	  bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);

	  std::string uniquePS = "Default";
	  std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	  if (chrName.substr(0,3) != "chr") chrName = std::string("chr").append(chrName);
	  int32_t regionChr = bam_name2id(hd, chrName.c_str());
	  int32_t regionStart = std::max(0, rec->pos - readLen);
	  int32_t regionEnd = (*svend) + readLen;
	  int32_t regionSize = regionEnd - regionStart;
	  typedef std::vector<uint32_t> TSplitPoints;
	  TSplitPoints sppH1;
	  sppH1.resize(regionSize, 0);
	  TSplitPoints sppH2;
	  sppH2.resize(regionSize, 0);
	  typedef std::vector<uint32_t> TBPCounts;
	  TBPCounts bpH1;
	  bpH1.resize(regionSize, 0);
	  TBPCounts bpH2;
	  bpH2.resize(regionSize, 0);

	  
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
	      } else uniquePS = rPS;

	      // Get the haplotype
	      uint8_t* hpptr = bam_aux_get(r, "HP");
	      if (hpptr) {
		int hap = bam_aux2i(hpptr);
		int32_t split = -1;
		if (splitPoint(r, split)) {
		  if ((split >= regionStart) && (split < regionEnd)) {
		    split -= regionStart;
		    if (hap==1) {
		      ++sppH1[split];
		      addBpCounts(r, regionStart, regionEnd, bpH1);
		    } else {
		      ++sppH2[split];
		      addBpCounts(r, regionStart, regionEnd, bpH2);
		    }
		  }
		}
	      }
	    }
	  }
	  bam_destroy1(r);
	  hts_itr_destroy(iter);

	  if (singlePhasedBlock) {
	    for(std::size_t i = 0; i<sppH1.size(); ++i) {
	      std::cout << sppH1[i] << ',' << sppH2[i] << "\t" << bpH1[i] << ',' << bpH2[i] << std::endl;
	    }
	  }
	  exit(-1);
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
