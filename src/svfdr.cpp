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
#include "htslib/vcf.h"


template<typename TVector>
inline void
_getMedian(TVector& v, typename TVector::value_type& med) {
  med = 0;
  if (v.size()) {
    typename TVector::iterator begin = v.begin();
    typename TVector::iterator end = v.end();
    std::nth_element(begin, begin + (end - begin) / 2, end);
    med = *(begin + (end - begin) / 2);
  }
}


int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " NA12878 <in.vcf.gz> <sample.bam>" << std::endl;
    return 1; 
  }

  std::string sampleName = std::string(argv[1]);

  // Load bcf file
  htsFile* ifile = bcf_open(argv[2], "r");
  if (!ifile) {
    std::cerr << "Fail to load " << argv[2] << "!" << std::endl;
    return 1;
  }
  
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  bcf1_t* rec = bcf_init();

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
	  std::cerr << bcf_hdr_id2name(hdr, rec->rid) << "\t" << (rec->pos + 1) << "\t" << *svend << "\t" << rec->d.id << std::endl;
	}
      }
    }
  }

  // Clean-up
  free(svend);
  free(svt);
  free(gt);
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  bcf_destroy(rec);

  return 0;
}
