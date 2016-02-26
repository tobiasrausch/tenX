/*
============================================================================
10X Phaser
============================================================================
Copyright (C) 2016 Tobias Rausch

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
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

struct Config {
  boost::filesystem::path outbam;
  boost::filesystem::path bamfile;
  boost::filesystem::path phasedbam;
  boost::filesystem::path variation;
};

struct Snp {
  uint32_t pos;
  char ref;
  char alt;
  uint32_t rch1;
  uint32_t rch2;
  uint32_t ach1;
  uint32_t ach2;
  uint32_t fc;
  Snp() : pos(0), ref('N'), alt('N'), rch1(0), rch2(0), ach1(0), ach2(0), fc(0) {}
  Snp(uint32_t p, char r, char a) : pos(p), ref(r), alt(a), rch1(0), rch2(0), ach1(0), ach2(0), fc(0) {}
};

template<typename TRecord>
struct SortSnps : public std::binary_function<TRecord, TRecord, bool> {
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    return s1.pos < s2.pos;
  }
};

inline uint32_t
lastAlignedPosition(bam1_t const* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  uint32_t alen = 0;
  for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL)) alen += bam_cigar_oplen(cigar[i]);
  return rec->core.pos + alen;
}

template<typename TSnpVector>
inline void
_refAltCount(bam1_t const* rec, TSnpVector& snps, int const hap) {
  typedef typename TSnpVector::value_type TSnp;
  
  // Annotate SNPs
  typename TSnpVector::iterator iSnp = std::lower_bound(snps.begin(), snps.end(), TSnp(rec->core.pos, 'A', 'A'), SortSnps<TSnp>());
  typename TSnpVector::iterator iSnpEnd = std::upper_bound(snps.begin(), snps.end(), TSnp(lastAlignedPosition(rec), 'A', 'A'), SortSnps<TSnp>());
  if (iSnp != iSnpEnd) {
    std::string sequence;
    sequence.resize(rec->core.l_qseq);
    uint8_t* seqptr = bam_get_seq(rec);
    for (std::size_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    uint32_t* cigar = bam_get_cigar(rec);
    for(;iSnp != iSnpEnd; ++iSnp) {
      int32_t gp = rec->core.pos; // Genomic position
      int32_t sp = 0; // Sequence position
      bool foundChar = false;
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	  if (gp + bam_cigar_oplen(cigar[i]) < iSnp->pos) {
	    gp += bam_cigar_oplen(cigar[i]);
	    sp += bam_cigar_oplen(cigar[i]);
	  } else {
	    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
	      if (gp == iSnp->pos) {
		foundChar = true;
		break;
	      }
	    }
	    if (foundChar) break;
	  }
	}
      }
      if (foundChar) {
	if (sequence[sp] == iSnp->ref) { 
	  if (hap == 1) ++iSnp->rch1;
	  else ++iSnp->rch2;
	} else if (sequence[sp] == iSnp->alt) {
	  if (hap == 1) ++iSnp->ach1;
	  else ++iSnp->ach2;
	} else ++iSnp->fc;
      }
    }
  }
}

template<typename TConfig, typename TGenomicSnps>
inline int32_t 
_loadMarkers(TConfig const& c, bam_hdr_t const* hdr, TGenomicSnps& snps) {
  typedef typename TGenomicSnps::value_type TSnpVector;
  typedef typename TSnpVector::value_type TSnp;

  // Open phased bam file
  samFile* phasedsam = sam_open(c.phasedbam.string().c_str(), "r");
  if (phasedsam == NULL) {
    std::cerr << "Fail to open file " << c.phasedbam.string() << std::endl;
    return 1;
  }
  hts_idx_t* phasedidx = sam_index_load(phasedsam, c.phasedbam.string().c_str());
  if (phasedidx == NULL) {
    std::cerr << "Fail to open index for " << c.phasedbam.string() << std::endl;
    return 1;
  }
  bam_hdr_t* phasedhdr = sam_hdr_read(phasedsam);
  if (phasedhdr == NULL) {
    std::cerr << "Fail to open header for " << c.phasedbam.string() << std::endl;
    return 1;
  }

  // Open VCF file
  htsFile* ifile = bcf_open(c.variation.string().c_str(), "r");
  if (ifile == NULL) {
    std::cerr << "SNP VCF files is missing " << c.variation.string() << std::endl;
    return 1;
  }
  hts_idx_t* bcfidx = bcf_index_load(c.variation.string().c_str());
  if (bcfidx == NULL) {
    std::cerr << "SNP VCF index file is missing " << c.variation.string() << std::endl;
    return 1;
  }
  bcf_hdr_t* vcfh = bcf_hdr_read(ifile);
    
  // Info
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parsing variation data" << std::endl;
  boost::progress_display show_progress(hdr->n_targets);

  // Load SNPs
  snps.clear();
  snps.resize(hdr->n_targets);
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    ++show_progress;
    // Collect Snps for this chromosome
    TSnpVector chrSnps;
    std::string chrName(hdr->target_name[refIndex]);
    uint32_t chrid = bcf_hdr_name2id(vcfh, chrName.c_str());
    hts_itr_t* itervcf = bcf_itr_queryi(bcfidx, chrid, 0, hdr->target_len[refIndex]);
    bcf1_t* var = bcf_init();
    while (bcf_itr_next(ifile, itervcf, var) >= 0) {
      bcf_unpack(var, BCF_UN_STR);
      std::vector<std::string> alleles;
      for(std::size_t i = 0; i<var->n_allele; ++i) alleles.push_back(std::string(var->d.allele[i]));
      // Only bi-allelic SNPs
      if ((alleles.size() == 2) && (alleles[0].size() == 1) && (alleles[1].size() == 1)) chrSnps.push_back(TSnp(var->pos, alleles[0][0], alleles[1][0]));
    }
    bcf_destroy(var);
    hts_itr_destroy(itervcf);
      
    // Sort Snps by position
    std::sort(chrSnps.begin(), chrSnps.end(), SortSnps<TSnp>());

    // Which Snps are informative?
    if (!chrSnps.empty()) {
      uint32_t tid = bam_name2id(phasedhdr, chrName.c_str());
      hts_itr_t* iter = sam_itr_queryi(phasedidx, tid, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(phasedsam, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	uint8_t *psptr = bam_aux_get(rec, "PS");
	if (psptr) {
	  // Only consider reads belonging to a phased block
	  char* ps = (char*) (psptr + 1);
	  std::string rPS = std::string(ps);
	  
	  // Get the haplotype
	  uint8_t* hpptr = bam_aux_get(rec, "HP");
	  if (hpptr) {
	    int hap = bam_aux2i(hpptr);
	    // Count
	    _refAltCount(rec, chrSnps, hap);
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    
      // Remove ambiguous SNPs
      for(typename TSnpVector::const_iterator itS = chrSnps.begin(); itS != chrSnps.end(); ++itS) {
	if (!itS->fc) {
	  if (((!itS->rch1) && (!itS->ach2)) || ((!itS->rch2) && (!itS->ach1))) {
	    if (((itS->rch1!=0) && (itS->ach2!=0)) || ((itS->rch2!=0) && (itS->ach1!=0))) {
	      //std::cerr << hdr->target_name[refIndex] << "\t" << itS->pos << "\t" << itS->ref << "\t" << itS->alt << "\t" << itS->rch1 << "\t" << itS->rch2 << "\t" << itS->ach1 << "\t" << itS->ach2 << "\t" << itS->fc << std::endl;
	      snps[refIndex].push_back(*itS);
	    }
	  }
	}
      }
      //std::cerr << hdr->target_name[refIndex] << "\t" << snps[refIndex].size() << std::endl;
    }
  }
  // Close VCF
  bcf_hdr_destroy(vcfh);
  hts_idx_destroy(bcfidx);
  bcf_close(ifile);
    
  return 0;
}



int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("phasebam.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("outbam,o", boost::program_options::value<boost::filesystem::path>(&c.outbam)->default_value("phased.bam"), "Phased output BAM file")
    ("phasedbam,p", boost::program_options::value<boost::filesystem::path>(&c.phasedbam), "Phased BAM file")
    ("variation,v", boost::program_options::value<boost::filesystem::path>(&c.variation), "SNP VCF file")
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
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -p <phased.bam> -v <snps.bcf> <unphased.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 

  // Check phased BAM file
  if (!(boost::filesystem::exists(c.phasedbam) && boost::filesystem::is_regular_file(c.phasedbam) && boost::filesystem::file_size(c.phasedbam))) {
    std::cerr << "Input phased BAM file is missing: " << c.phasedbam.string() << std::endl;
    return 1;
  }

  // Check input BAM file
  if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
    std::cerr << "Input phased BAM file is missing: " << c.bamfile.string() << std::endl;
    return 1;
  }
  
  // Check variation VCF file
  if (!(boost::filesystem::exists(c.variation) && boost::filesystem::is_regular_file(c.variation) && boost::filesystem::file_size(c.variation))) {
    std::cerr << "Input SNP VCF file is missing: " << c.variation.string() << std::endl;
    return 1;
  }


  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Load bam files
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  if (samfile == NULL) {
    std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
    return 1;
  }
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  if (idx == NULL) {
    std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
    return 1;
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile);
  if (hdr == NULL) {
    std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
    return 1;
  }

  // Load variation data
  typedef std::vector<Snp> TSnpVector;
  typedef std::vector<TSnpVector> TGenomicSnps;
  TGenomicSnps snps;
  int32_t r = _loadMarkers(c, hdr, snps);
  if (r) return r;

  // Open output file
  samFile* outbam = sam_open(c.outbam.string().c_str(), "wb");
  if (outbam == NULL) {
    std::cerr << "Fail to open file " << c.outbam.string() << std::endl;
    return 1;
  }
  sam_hdr_write(outbam, hdr);

  // Assign reads to SNPs
  uint32_t assignedReads = 0;
  uint32_t unassignedReads = 0;
  uint32_t assignedBases = 0;
  uint32_t unassignedBases = 0;
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Phasing reads" << std::endl;
  boost::progress_display show_progress(hdr->n_targets);
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    ++show_progress;
    if (snps[refIndex].empty()) continue;

    hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(samfile, iter, rec) >= 0) {

      uint32_t hp1votes = 0;
      uint32_t hp2votes = 0;
      TSnpVector::const_iterator iSnp = std::lower_bound(snps[refIndex].begin(), snps[refIndex].end(), Snp(rec->core.pos, 'A', 'A'), SortSnps<Snp>());
      TSnpVector::const_iterator iSnpEnd = std::upper_bound(snps[refIndex].begin(), snps[refIndex].end(), Snp(lastAlignedPosition(rec), 'A', 'A'), SortSnps<Snp>());
      if (iSnp != iSnpEnd) {
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	uint32_t* cigar = bam_get_cigar(rec);
	for(;iSnp != iSnpEnd; ++iSnp) {
	  uint32_t gp = rec->core.pos; // Genomic position
	  uint32_t sp = 0; // Sequence position
	  bool foundChar = false;
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	      if (gp + bam_cigar_oplen(cigar[i]) < iSnp->pos) {
		gp += bam_cigar_oplen(cigar[i]);
		sp += bam_cigar_oplen(cigar[i]);
	      } else {
		for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		  if (gp == iSnp->pos) {
		    foundChar = true;
		    break;
		  }
		}
		if (foundChar) break;
	      }
	    }
	  }
	  if (foundChar) {
	    if (sequence[sp] == iSnp->ref) { 
	      if (iSnp->rch1) ++hp1votes;
	      else ++hp2votes;
	    } else if (sequence[sp] == iSnp->alt) {
	      if (iSnp->ach1) ++hp1votes;
	      else ++hp2votes;
	    }
	  }
	}
      }
      int32_t hp = 0;
      if (hp1votes > 2*hp2votes) hp = 1;
      else if (hp2votes > 2*hp1votes) hp = 2;
      if (hp) {
	++assignedReads;
	assignedBases += rec->core.l_qseq;
	bam_aux_append(rec, "HP", 'i', 4, (uint8_t*)&hp);
	sam_write1(outbam, hdr, rec);
      } else {
	++unassignedReads;
	unassignedBases += rec->core.l_qseq;
      }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
  }

  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  // Statistics
  uint32_t sumReads = assignedReads + unassignedReads;
  uint32_t sumBases = assignedBases + unassignedBases;
  std::cout << "Assigned Reads=" << assignedReads << ", Unassigned Reads=" << unassignedReads << ", Fraction assigned=" << (float) assignedReads / (float) sumReads << std::endl;
  std::cout << "Assigned Bases=" << assignedBases << ", Unassigned Bases=" << unassignedBases << ", Fraction assigned=" << (float) assignedBases / (float) sumBases << std::endl;


#ifdef PROFILE
  ProfilerStop();
#endif


  return 0;
}
