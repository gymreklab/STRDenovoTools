/*
This file is modified from Thomas Willems's HipSTR: https://github.com/tfwillems/HipSTR
*/


#include <set>
#include <fstream>
#include <math.h>
#include <sys/stat.h>

#include "vcf_reader.h"

std::string GRID_INFO_TAG         = "GRID";

namespace VCF {

  const std::vector<std::string>& Variant::get_samples() const {
    return vcf_reader_->get_samples();
  }

  int Variant::round_allele_length(int allele) const {
    int period = 0;
    get_INFO_value_single_int("PERIOD", period);
    int rounded_allele = round(static_cast<float>(allele)/static_cast<float>(period))*period;
    return rounded_allele;
  }

  int Variant::num_alleles_by_length(const bool& round_allele) const {
    int ref_allele_size = (int)alleles_.front().size();
    std::set<int> lengths;
    for (auto iter=alleles_.begin(); iter!=alleles_.end(); iter++) {
      if (round_allele) {
	lengths.insert(round_allele_length((int)iter->size()-ref_allele_size));
      } else {
	lengths.insert((int)iter->size());
      }
    }
    return (int)lengths.size();
  }

  float Variant::heterozygosity() const {
    std::vector<double> allele_freqs(num_alleles(), 1.0); // Use a one sample pseudocount
    double total_count = num_alleles();
    // Get counts
    const std::vector<std::string> samples = get_samples();
    int gt_a, gt_b;
    for (auto iter = samples.begin(); iter != samples.end(); iter++) {
      if (sample_call_missing(*iter)) {
	continue;
      }
      get_genotype(*iter, gt_a, gt_b);
      allele_freqs[gt_a]++;
      allele_freqs[gt_b]++;
      total_count += 2;
    }
    // Compute heterozygosity
    double x = 0.0;
    for (int i = 0; i < allele_freqs.size(); i++) {
      double p = allele_freqs[i]/total_count;
      x += p*p;
    }
    return 1-x; 
  }

  float Variant::heterozygosity_by_length(const bool& round_alleles) const {
    std::vector<double> allele_freqs(num_alleles_by_length(round_alleles), 1.0);
    double total_count = num_alleles_by_length(round_alleles);
    const std::vector<std::string> samples = get_samples();
    int gt_a, gt_b;
    for (auto iter = samples.begin(); iter != samples.end(); iter++) {
      if (sample_call_missing(*iter)) {
	continue;
      }
      get_genotype(*iter, gt_a, gt_b);
      allele_freqs[GetLengthIndexFromGT(gt_a)]++;
      allele_freqs[GetLengthIndexFromGT(gt_b)]++;
      total_count += 2;
    }
    // Compute heterozygosity
    double x = 0.0;
    for (int i = 0; i < allele_freqs.size(); i++) {
      double p = allele_freqs[i]/total_count;
      x += p*p;
    }
    return 1-x; 
  }

  int Variant::num_gangstr_alleles() const {
    std::vector<std::int32_t> grid_vals;
    if (!has_info_field(GRID_INFO_TAG)) {
      PrintMessageDieOnError("Required INFO field " + GRID_INFO_TAG + " not present in VCF", M_ERROR);
    }
    get_INFO_value_multiple_ints(GRID_INFO_TAG, grid_vals);
    if (grid_vals.size() != 2) {
      PrintMessageDieOnError("Inproperly formatted GRID field. Should have 2 values", M_ERROR);
    }
    return grid_vals[1]-grid_vals[0]+1;
  }

  void Variant::build_alleles_by_length(const bool& round_alleles) {
    //    std::cerr << "building allele length map " << std::endl;
    // Build map of GT->allele length. Assume alleles ordered by length
    // Except reference allele, which is always first
    const std::vector<std::string> alleles = get_alleles();
    int ref_allele_size = (int) alleles.front().size();
    length_allele_sizes_.resize(num_alleles_by_length(round_alleles), 0);
    // First element is always ref allele
    length_allele_sizes_[0] = 0;
    gt_to_length_index_[0] = 0;
    int max_index = 0;
    // Go through rest of alleles.
    int allele_index = 1;
    int allele_size;
    // First check the first one
    if (alleles.size() > 1) {
      allele_size = (int)alleles[1].size() - ref_allele_size;
      if (round_alleles) {
	allele_size = round_allele_length(allele_size);
      }
      if (allele_size == 0) {
	allele_index = 0;
      }
    }
    int prev_index = 0; // Only gets set when we encounter ref allele
    for (int i = 1; i < alleles.size(); i++) {
      int len = (int)alleles[i].size() - ref_allele_size;
      if (round_alleles) {
	len = round_allele_length(len);
      }
      assert(len >= allele_size); // Assume alleles ordered by length
      if (len == 0) { // If this alleles is ref length, keep track of prev index
	allele_size = len;
	if (allele_index != 0) {
	  prev_index = allele_index;
	}
	allele_index = 0;
      } else if (allele_size == 0) { // If previous allele was ref length, continue indexing where we left off
	allele_size = len;
	allele_index = prev_index + 1;
      } else if (len == allele_size) {
	allele_size = len;
      } else {
	allele_size = len;
	allele_index++;
      }
      gt_to_length_index_[i] = allele_index;
      length_allele_sizes_[allele_index] = allele_size;
      if (allele_index > max_index) {
	max_index = allele_index;
      }
    }
    assert(max_index == num_alleles_by_length(round_alleles)-1);
  }

  int Variant::GetLengthIndexFromGT(const int& gt_index) const {
    return gt_to_length_index_.at(gt_index);
  }

  int Variant::GetSizeFromLengthAllele(const int& length_gt_index) const {
    return length_allele_sizes_[length_gt_index];
  }

  void Variant::get_length_genotype(int sample_index, int& gt_al, int& gt_bl) const {
    int gt_a, gt_b; // First get regular genotype
    get_genotype(sample_index, gt_a, gt_b);
    // Now convert to lengths (relative to reference)
    int lenind_a = GetLengthIndexFromGT(gt_a);
    int lenind_b = GetLengthIndexFromGT(gt_b);
    gt_al = GetSizeFromLengthAllele(lenind_a);
    gt_bl = GetSizeFromLengthAllele(lenind_b);
  }

  void Variant::get_length_genotype(const std::string& sample, int& gt_al, int& gt_bl) const {
    int gt_a, gt_b; // First get regular genotype
    get_genotype(sample, gt_a, gt_b);
    // Now convert to lengths (relative to reference)
    int lenind_a = GetLengthIndexFromGT(gt_a);
    int lenind_b = GetLengthIndexFromGT(gt_b);
    gt_al = GetSizeFromLengthAllele(lenind_a);
    gt_bl = GetSizeFromLengthAllele(lenind_b);
  }

  void Variant::get_genotype(const std::string& sample, int& gt_a, int& gt_b) const {
    int sample_index = vcf_reader_->get_sample_index(sample);
    if (sample_index == -1)
      gt_a = gt_b = -1;
    else {
      gt_a = gt_1_[sample_index];
      gt_b = gt_2_[sample_index];
    }
  }

  bool Variant::sample_call_missing(const std::string& sample) const {
    int sample_index = vcf_reader_->get_sample_index(sample);
    return (sample_index == -1 ? true : missing_[sample_index]);
  }

  void Variant::extract_alleles(){
    for (int i = 0; i < vcf_record_->n_allele; i++)
      alleles_.push_back(vcf_record_->d.allele[i]);
  }

  void Variant::extract_genotypes(){
    int   mem = 0;
    int* gts_ = NULL;
    std::string GT_KEY = "GT";
    int num_entries = bcf_get_format_int32(vcf_header_, vcf_record_, GT_KEY.c_str(), &gts_, &mem);
    if (num_entries <= 0)
      PrintMessageDieOnError("Failed to extract the genotypes from the VCF record", M_ERROR);
    if (num_entries != num_samples_ && num_entries != 2*num_samples_)
      PrintMessageDieOnError("Incorrect number of genotypes extracted from the VCF record", M_ERROR);

    if (num_entries == num_samples_){
      // When all sample genotypes are missing, the number of extracted GT's is the
      // same as the number of samples
      missing_ = std::vector<bool>(num_samples_, true);
      phased_  = std::vector<bool>(num_samples_, false);
      gt_1_    = std::vector<int>(num_samples_, -1);
      gt_2_    = std::vector<int>(num_samples_, -1);
    }
    else {
      missing_.reserve(num_samples_);
      phased_.reserve(num_samples_);
      gt_1_.reserve(num_samples_);
      gt_2_.reserve(num_samples_);
    
      int gt_index = 0;
      for (int i = 0; i < num_samples_; i++){
	if (bcf_gt_is_missing(gts_[gt_index]) || bcf_gt_is_missing(gts_[gt_index+1])){
	  missing_.push_back(true);
	  phased_.push_back(false);
	  gt_1_.push_back(-1);
	  gt_2_.push_back(-1);
	}
	else {
	  missing_.push_back(false);
	  phased_.push_back(bcf_gt_is_phased(gts_[gt_index+1]));
	  gt_1_.push_back(bcf_gt_allele(gts_[gt_index]));
	  // For haploid, set two alleles to be the same.
	  // This check based on https://github.com/samtools/htslib/blob/cb8a05138e3904f7993d71026d3ee463e0cc80e9/htslib/vcf.h#L1296
	  // TODO check this. kind of sketchy bc numbers don't exactly match
	  // bcf_int32_vector_end (-2147483647). 
	  if ((bcf_gt_allele(gts_[gt_index+1])+1)*2+1 != bcf_int32_vector_end) {
	    // diploid
	    gt_2_.push_back(bcf_gt_allele(gts_[gt_index+1]));
	  } else {
	    // haploid
	    gt_2_.push_back(bcf_gt_allele(gts_[gt_index]));
	  }
	}
	gt_index += 2;
      }
    }
    free(gts_);
  }

  void VCFReader::open(const std::string& filename){
    const char* cfilename = filename.c_str();
    
    std::ifstream checkfile(cfilename);
    if (!checkfile.good()) {
      PrintMessageDieOnError("VCF file could not be opened " + filename, M_ERROR);
    }

    if (bgzf_is_bgzf(cfilename) != 1)
      PrintMessageDieOnError("VCF file is not in a valid bgzipped file. Please ensure that bgzip was used to compress it", M_ERROR);
  
    char *fnidx = (char*) calloc(strlen(cfilename) + 5, 1);
    strcat(strcpy(fnidx, cfilename), ".tbi");
    struct stat stat_tbi, stat_vcf;
    stat(fnidx, &stat_tbi);
    stat(cfilename, &stat_vcf);
    if (stat_vcf.st_mtime > stat_tbi.st_mtime)
      PrintMessageDieOnError("The tabix index for the VCF file is older than the VCF itself. Please reindex the VCF with tabix", M_ERROR);
    free(fnidx);
  
    if ((vcf_input_ = hts_open(cfilename, "r")) == NULL)
      PrintMessageDieOnError("Failed to open the VCF file", M_ERROR);
    if ((tbx_input_ = tbx_index_load(cfilename)) == NULL)
      PrintMessageDieOnError("Failed to open the VCF file's tabix index", M_ERROR);
  
    int nseq;
    const char** seq = tbx_seqnames(tbx_input_, &nseq);
    for (int i = 0; i < nseq; i++)
      chroms_.push_back(seq[i]);
    free(seq);
  
    if (chroms_.size() == 0)
      PrintMessageDieOnError("VCF does not contain any chromosomes", M_ERROR);
  
    vcf_header_  = bcf_hdr_read(vcf_input_);
    tbx_iter_    = tbx_itr_querys(tbx_input_, chroms_.front().c_str());
    chrom_index_ = 0;
  
    for (int i = 0; i < bcf_hdr_nsamples(vcf_header_); i++){
      samples_.push_back(vcf_header_->samples[i]);
      sample_indices_[vcf_header_->samples[i]] = i;
    }
  }

  bool VCFReader::get_next_variant(Variant* variant, const bool& round){
    vcf_record_ = bcf_init(); // need to destroy after variant goes out of scope
    if ((tbx_iter_ != NULL) && tbx_itr_next(vcf_input_, tbx_input_, tbx_iter_, &vcf_line_) >= 0){
      if (vcf_parse(&vcf_line_, vcf_header_, vcf_record_) < 0)
	PrintMessageDieOnError("Failed to parse VCF record", M_ERROR);
      *variant = Variant(vcf_header_, vcf_record_, this, round);
      return true;
    }
  
    if (jumped_) {
      bcf_destroy(vcf_record_);
      return false;
    }
  
    while (chrom_index_+1 < chroms_.size()){
      chrom_index_++;
      tbx_itr_destroy(tbx_iter_);
      tbx_iter_ = tbx_itr_querys(tbx_input_, chroms_[chrom_index_].c_str());
    
      if ((tbx_iter_ != NULL) && tbx_itr_next(vcf_input_, tbx_input_, tbx_iter_, &vcf_line_) >= 0){
	if (vcf_parse(&vcf_line_, vcf_header_, vcf_record_) < 0)
	  PrintMessageDieOnError("Failed to parse VCF record", M_ERROR);
	*variant = Variant(vcf_header_, vcf_record_, this, round);
	return true;
      }
    }
    bcf_destroy(vcf_record_);
    return false;
  }

}
