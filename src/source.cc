#define VERSION "v0.9.0"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include <fstream>
#include <limits>
#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

#include <plinkio/plinkio.h>

#define DEFAULT_SIGSQ 1.0f
#define DEFAULT_CAUSAL_PI 0.001f
#define DEFAULT_HSQ 0.7f
#define DEFAULT_RG 0.0f
#define DEFAULT_K 0.1f

struct SimuOptions {
  std::string bfile;
  std::string bfile_chr;
  bool qt;
  bool cc;
  int num_traits;
  std::vector<float> k;
  std::vector<int> ncas;
  std::vector<int> ncon;
  std::vector<float> hsq;
  int num_components;
  std::vector<std::string> causal_variants;
  std::vector<int> causal_n;
  std::vector<float> causal_pi;
  std::vector<std::string> causal_regions;
  std::vector<float> trait1_sigsq;
  std::vector<float> trait2_sigsq;
  std::vector<float> rg;
  bool gcta_sigma;
  
  // see documentation at github about --norm-effect flag;
  // however, beware that everywhere in this file, regardless of this flag, the computations are done w.r.t. 0-1-2 genotypes (additively coded).
  // it implies that if user provides effect sizes, and --norm-effect is on, we first convert user effect sizes into effect sizes w.r.t. 0-1-2, and convert them back into per-normalized-genotypes effect sizes before saving the resulting .causals file 
  bool norm_effect;

  std::string out;
  boost::int64_t seed;
  bool verbose;

  // derived options
  std::vector<std::string> bfiles;
  int num_samples;
  int num_variants;
  std::string out_pheno;
  std::vector<std::string> out_causals;

  SimuOptions() {
    num_samples = -1;
    num_variants = 0;
    boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
    seed = (boost::posix_time::microsec_clock::local_time() - time_epoch).ticks();
    verbose = false;
  }
};

class Logger {
 public:
  Logger(std::string path) : log_file_(path) { 
    // boost::posix_time::time_facet *facet = new boost::posix_time::time_facet("%d-%b-%Y %H:%M:%S");
    // std::cout.imbue(std::locale(std::cout.getloc(), facet));
    // log_file_.imbue(std::locale(log_file_.getloc(), facet));
    // std::cout << std::setprecision(9);
  }

  template <typename T>
  Logger& operator<< ( const T& rhs) {
    std::cout << rhs;
    log_file_ << rhs;
    return *this;
  }

  template <typename T>
  Logger& operator<< (const std::vector<T>& rhs) {
    for (int i = 0; i < rhs.size(); i++)
      (*this) << ((i>0) ? " " : "") << rhs[i];
    return (*this);
  }

 private:
  std::ofstream log_file_;
};

void log_header(int argc, char *argv[], Logger &log) {
  std::string header(
"*********************************************************************\n"
"* SIMU - library for simulation of GWAS summary statistics\n"
"* Version " VERSION "\n"
"* (C) 2018 Oleksandr Frei\n"
"* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
"* GNU General Public License v3\n"
"*********************************************************************\n");
 
  log << header;
  if (argc==0) return;
  log  << "Call:\n" << argv[0] << " ";
  for (int i = 1; i < argc; i++) {
    if (strlen(argv[i]) == 0) continue;
    if (argv[i][0] == '-') log << "\\\n\t";
    log << argv[i] << " ";
  }
  log << "\n\n";
}

// https://stackoverflow.com/questions/4107087/accepting-negative-doubles-with-boostprogram-options
std::vector<po::option> ignore_numbers(std::vector<std::string>& args) {
  std::vector<po::option> result;
  int pos = 0;
  while(!args.empty()) {
    const auto& arg = args[0];

    bool is_number = false;
    try {
      boost::lexical_cast<double>(arg);
      is_number = true;
    } catch(...) {
    }

    if(is_number) {
      result.push_back(po::option());
      po::option& opt = result.back();

      opt.position_key = pos++;
      opt.value.push_back(arg);
      opt.original_tokens.push_back(arg);

      args.erase(args.begin());
    } else {
      break;
    }
  }

  return result;
}

class PioFile {
 public:
  PioFile(std::string bfile) : bfile_(bfile) {
    if( pio_open( &plink_file_, bfile.c_str() ) != PIO_OK )
    {
      std::stringstream ss; ss << "ERROR: Could not open " << bfile;
      throw std::runtime_error(ss.str());
    }

    if( !pio_one_locus_per_row( &plink_file_ ) )
    {
      std::stringstream ss; ss << "ERROR: Unsupported format in " << bfile << ", snps must be rows, samples must be columns";
      throw std::runtime_error(ss.str());
    }

    row_size_ = pio_row_size( &plink_file_ );
    if (row_size_ == 0) {
      std::stringstream ss; ss << "ERROR: Unsupported format in " << bfile << ", rows appears to have zero size";
      throw std::runtime_error(ss.str());
    }

    num_samples_ = pio_num_samples( &plink_file_ );
    num_loci_ = pio_num_loci( &plink_file_ );
    cur_locus_ = 0;
  }

  pio_status_t skip_row() { 
    pio_status_t retval = pio_skip_row(&plink_file_);
    if (retval == PIO_ERROR) {
      std::stringstream ss; ss << "ERROR: Error reading from " << bfile_;
      throw std::runtime_error(ss.str());
    }
    cur_locus_++;
    return retval;
  }

  pio_status_t next_row(std::vector<snp_t> *buffer) {
    if (buffer->size() != row_size_) buffer->resize(row_size_);
    pio_status_t retval = pio_next_row(&plink_file_, &(*buffer)[0]);
    if (retval == PIO_ERROR) {
      std::stringstream ss; ss << "ERROR: Error reading from " << bfile_;
      throw std::runtime_error(ss.str());
    }
    cur_locus_++;
    return retval;
  }

  void reset_row() {
    pio_reset_row(&plink_file_);
    cur_locus_ = 0;
  }

  pio_sample_t* get_sample(size_t sample_id) {
    return pio_get_sample(&plink_file_, sample_id);
  }

  pio_locus_t* get_locus(size_t locus_id) {
    return pio_get_locus(&plink_file_, locus_id);
  }

  // return loci that will be read by the next call to next_row
  pio_locus_t* get_cur_locus() {
    return pio_get_locus(&plink_file_, cur_locus_);
  }

  int row_size() const { return row_size_; }
  int num_samples() const { return num_samples_; }
  int num_loci() const { return num_loci_; }  // I use "loci" here to be consistent with plinkio terminology, but everywhere else it goes as "variants"

  ~PioFile() {
    pio_close( &plink_file_ );
  }

 private:
  std::string bfile_;
  pio_file_t plink_file_;
  int row_size_;
  int num_samples_;
  int num_loci_;
  int cur_locus_;
};

class PioFiles {
 public:
  PioFiles() { cur_file_ = 0; num_variants_ = 0; }
  void push_back(std::shared_ptr<PioFile> pio_file) {
    pio_files_.push_back(pio_file);
    num_variants_ += pio_file->num_loci();
  }

  pio_status_t skip_row() {
    while (true) {
      if (cur_file_ >= pio_files_.size()) {
        return PIO_END;
      }

      if (pio_files_[cur_file_]->skip_row() == PIO_OK)
        return PIO_OK;
       
      cur_file_++;
    }
  }

  pio_status_t next_row(std::vector<snp_t>* buffer) {
    while (true) {
      if (cur_file_ >= pio_files_.size()) {
        return PIO_END;
      }

      if (pio_files_[cur_file_]->next_row(buffer) == PIO_OK)
        return PIO_OK;
       
      cur_file_++;
    }
  }

  void reset_row() {
    for (int i = 0; i < pio_files_.size(); i++)
      pio_files_[i]->reset_row();
    cur_file_ = 0;
  }

  pio_sample_t* get_sample(size_t sample_id) {
    return pio_files_[0]->get_sample(sample_id);
  }

  pio_locus_t* get_locus(size_t locus_id) {
    for (int i = 0; i < pio_files_.size(); i++) {
      if (locus_id < pio_files_[i]->num_loci())
        return pio_files_[i]->get_locus(locus_id);
      locus_id -= pio_files_[i]->num_loci();
    }
    return nullptr;
  }

  // return loci that will be read by the next call to next_row
  pio_locus_t* get_cur_locus() {
    return pio_files_[cur_file_]->get_cur_locus();
  }

  void find_variant_index_map() {
    if (!snp_to_id_.empty())
      throw(std::runtime_error("internal error - find_variant_index_map must be called once"));

    for (int variant_index = 0; variant_index < num_variants_; variant_index++) {
      pio_locus_t* locus = this->get_locus(variant_index);
      if (locus == nullptr) {
        std::stringstream ss; ss << "ERROR: unable to read locus with index " << variant_index;
        throw(std::runtime_error(ss.str()));
      }
      auto retval = snp_to_id_.insert({locus->name, variant_index});
      if (retval.second == false) {
        std::stringstream ss; ss << "ERROR: duplicated variant name: " << locus->name;
        throw(std::runtime_error(ss.str()));
      }
    }
  }

  int get_locus_id(std::string locus_name) {
    auto it = snp_to_id_.find(locus_name);
    if (it == snp_to_id_.end()) {
      std::stringstream ss; ss << "ERROR: locus with name " << locus_name << " does not exist in the input; check your --bfile and/or --causal-variants, --causal-regions arguments";
      throw(std::runtime_error(ss.str()));
    }
    return it->second;
  }

  int num_variants() const { return num_variants_; }

 private:
  std::vector<std::shared_ptr<PioFile>> pio_files_;
  std::map<std::string, int> snp_to_id_;
  int cur_file_;
  int num_variants_;
};

void read_causal_variants_file(std::string file, Logger* log,
                               PioFiles* pio_files,
                               std::vector<int>* variant_indices,
                               std::vector<double>* effect1_per_variant,
                               std::vector<double>* effect2_per_variant) {
  std::ifstream infile(file);
  std::string line;
  variant_indices->clear(); 
  if (effect1_per_variant != nullptr) effect1_per_variant->clear();
  if (effect2_per_variant != nullptr) effect2_per_variant->clear();
  int num_tokens = -1;
  int line_num = 1;
  if (log != nullptr) (*log) << "Reading '" << file << "'... ";
  while (std::getline(infile, line))
  {
    std::vector<std::string> tokens;
    const std::string separators = " \t\n\r";
    boost::trim_if(line, boost::is_any_of(separators));
    boost::split(tokens, line, boost::is_any_of(separators), boost::token_compress_on);
    if (line.empty() || tokens.size() == 0) continue;
    if (num_tokens == -1) num_tokens = tokens.size();
    if (tokens.size() != num_tokens) {
      std::stringstream ss; ss << "ERROR: failed to parse line " << line_num << " from '" << file << "', invalid number of tokens (" << tokens.size() << " found, " << num_tokens << " expected)";
      throw std::runtime_error(ss.str()); 
    }

    int variant_index = pio_files->get_locus_id(tokens[0]);
    if ((variant_index == -1) || (variant_index > pio_files->num_variants())) {
      std::stringstream ss; ss << "ERROR: can not find variant '" << tokens[0] << "' in input --bfile / --bfile-chr";
      throw std::runtime_error(ss.str()); 
    }

    double effect1, effect2;
    try {
      if (tokens.size() > 1) effect1 = boost::lexical_cast<double>(tokens[1]);
    } catch(...) {
      std::stringstream ss; ss << "ERROR: failed to parse line " << line_num << " from '" << file << "', effect size '" << tokens[1] << "' is not a number";
      throw std::runtime_error(ss.str());
    }

    try {
      if (tokens.size() > 2) effect2 = boost::lexical_cast<double>(tokens[2]);
    } catch(...) {
      std::stringstream ss; ss << "ERROR: failed to parse line " << line_num << " from '" << file << "', effect size '" << tokens[2] << "' is not a number";
      throw std::runtime_error(ss.str());
    }

    variant_indices->push_back(variant_index);
    if ((effect1_per_variant != nullptr) && (tokens.size() > 1)) effect1_per_variant->push_back(effect1);
    if ((effect2_per_variant != nullptr) && (tokens.size() > 2)) effect2_per_variant->push_back(effect2);
    line_num++;
  }

  if (variant_indices->empty()) {
    std::stringstream ss; ss << "ERROR: '" << file << "' is empty or contains no variants";
    throw std::runtime_error(ss.str());
  }
  if (log != nullptr) (*log) << variant_indices->size() << " variants found.\n";
}

void fix_and_validate(SimuOptions& simu_options, po::variables_map& vm, Logger& log,
                      PioFiles* pio_files) {
  // Validate and pre-process --bfile / --bfile-chr options.
  if (simu_options.bfile.empty() && simu_options.bfile_chr.empty())
    throw std::invalid_argument(std::string("ERROR: Either --bfile or --bfile-chr must be specified"));
  if (!simu_options.bfile_chr.empty()) {
    if (!boost::contains(simu_options.bfile_chr, "@"))
      simu_options.bfile_chr.append("@");
    for (int chri = 1; chri <= 22; chri++) {
      simu_options.bfiles.push_back(simu_options.bfile_chr);
      boost::replace_all(simu_options.bfiles.back(), "@", boost::lexical_cast<std::string>(chri));
    }
  } else {
    simu_options.bfiles.push_back(simu_options.bfile);
  }

  // Validate --cc and --qt options
  if (!simu_options.cc && !simu_options.qt)
    throw std::invalid_argument(std::string("ERROR: Either --qt or --cc must be specified"));

  // validate number of traits
  if (simu_options.num_traits <= 0 || simu_options.num_traits >= 3)
    throw std::invalid_argument(std::string("ERROR: --num-traits must be either 1 or 2"));

  // validate number of components
  if (simu_options.num_components <= 0)
    throw std::invalid_argument(std::string("ERROR: --num-components must be non-negative number"));

  // Validate --out options
  if (simu_options.out.empty())
    throw std::invalid_argument(std::string("ERROR: --out option must be specified"));
  simu_options.out_pheno = simu_options.out + ".pheno";
  for (int i = 0; i < simu_options.num_traits; i++)
    simu_options.out_causals.push_back(simu_options.out + "." + boost::lexical_cast<std::string>(i+1) + ".causals");

  if ( boost::filesystem::exists( simu_options.out_pheno ) )
    log << "WARNING: Target file " << simu_options.out_pheno << " already exists and will be overwritten\n"; 

  if (simu_options.qt && (vm.count("k") > 0 || vm.count("ncas") > 0 || vm.count("ncon") > 0)) {
    log << "WARNING: Options --k, --ncas, --ncon are not relevant to --qt, and will be ignored\n";
    simu_options.k.clear();
    simu_options.ncas.clear();
    simu_options.ncon.clear();
  }
  if ((simu_options.num_traits==1) && (vm.count("trait2-sigsq") > 0 || (vm.count("rg") > 0))) {
    log << "WARNING: Options --trait2-sigsq and --rg are not relevant for a single-trait simulations, and will be ignored\n";
    simu_options.trait2_sigsq.clear();
    simu_options.rg.clear();
  }

  if (simu_options.cc && (simu_options.ncas.empty() != simu_options.ncon.empty()))
    throw std::invalid_argument("ERROR: Options --ncas and --ncon must be either both present, or both absent");

  if (simu_options.cc) {
    if (!simu_options.k.empty() && (simu_options.k.size() != simu_options.num_traits))
      throw std::invalid_argument(std::string("ERROR: Number of --k values does not match the number of traits"));
    if (!simu_options.ncas.empty() && (simu_options.ncas.size() != simu_options.num_traits))
      throw std::invalid_argument(std::string("ERROR: Number of --ncas values does not match the number of traits"));
    if (!simu_options.ncon.empty() && (simu_options.ncon.size() != simu_options.num_traits))
      throw std::invalid_argument(std::string("ERROR: Number of --ncon values does not match the number of traits"));
  }
  if (!simu_options.hsq.empty() && (simu_options.hsq.size() != simu_options.num_traits))
    throw std::invalid_argument(std::string("ERROR: Number of --hsq values does not match the number of traits"));

  if (simu_options.cc) {
    // initialize default value for k (prevalence for case/control trait) and check for out of range values
    if (simu_options.k.empty())
      for (int i = 0; i < simu_options.num_traits; i++) simu_options.k.push_back(DEFAULT_K);
    for (auto val: simu_options.k) if (val <= 0 || val >= 1)
      throw std::invalid_argument(std::string("ERROR: Prevalence --k must be between 0.0 and 1.0"));
  }

  // inidialize default value for hsq (heritability) and check for out of range values
  if (simu_options.hsq.empty())
    for (int i = 0; i < simu_options.num_traits; i++) simu_options.hsq.push_back(DEFAULT_HSQ);
  for (auto val: simu_options.hsq) if (val < 0 || val > 1)
    throw std::invalid_argument(std::string("ERROR: Heritability --hsq must be between 0.0 and 1.0"));

  // initialize default value for --causal-pi, --causal-n, --causal-variants and check for out of range values
  if (!simu_options.causal_pi.empty() && (simu_options.causal_pi.size() != simu_options.num_components))
    throw std::invalid_argument(std::string("ERROR: Number of --causal-pi values does not match the number of components"));
  if (!simu_options.causal_n.empty() && (simu_options.causal_n.size() != simu_options.num_components))
    throw std::invalid_argument(std::string("ERROR: Number of --causal-n values does not match the number of components"));
  if (!simu_options.causal_variants.empty() && (simu_options.causal_variants.size() != simu_options.num_components))
    throw std::invalid_argument(std::string("ERROR: Number of --causal-variants values does not match the number of components"));
  if (((int)simu_options.causal_pi.empty() + (int)simu_options.causal_n.empty() + (int)simu_options.causal_variants.empty()) < 2)
    throw std::invalid_argument(std::string("ERROR: At most one of --causal-pi, --causal-n, --causal-variants options must be used"));
  if (simu_options.causal_pi.empty() && simu_options.causal_n.empty() && simu_options.causal_variants.empty())
    for (int i = 0; i < simu_options.num_components; i++) simu_options.causal_pi.push_back(DEFAULT_CAUSAL_PI);
  for (auto val: simu_options.causal_pi) if (val <= 0 || val > 1)
    throw std::invalid_argument(std::string("ERROR: --causal-pi values must be between 0.0 and 1.0"));
  for (auto val: simu_options.causal_n) if (val <= 0)
    throw std::invalid_argument(std::string("ERROR: --causal-n values must be a positive"));

  if (!simu_options.causal_regions.empty() && !simu_options.causal_variants.empty()) {
    log << "WARNING: Option --causal-regions can not be used together with --causal-variants, and will be ignored\n";
    simu_options.causal_regions.clear();
  }  
  if (!simu_options.causal_regions.empty() && (simu_options.causal_regions.size() != simu_options.num_components))
    throw std::invalid_argument(std::string("ERROR: Number of --causal-regions values does not match the number of components"));

  for (auto val: simu_options.causal_regions) if (!boost::filesystem::exists( val )) {
    std::stringstream ss; ss << "ERROR: input file " << val << " does not exist";
    throw std::runtime_error(ss.str());
  }

  for (auto val: simu_options.causal_variants) if (!boost::filesystem::exists( val )) {
    std::stringstream ss; ss << "ERROR: input file " << val << " does not exist";
    throw std::runtime_error(ss.str());
  } 

  // initialize default value for --trait1-sigsq and check for out of range values
  {
    auto& target = simu_options.trait1_sigsq;
    if (!target.empty() && (target.size() != simu_options.num_components))
      throw std::invalid_argument(std::string("ERROR: Number of --trait1-sigsq values does not match the number of components"));
    if (target.empty())
      for (int i = 0; i < simu_options.num_components; i++) target.push_back(DEFAULT_SIGSQ);
    for (auto val: target) if (val < 0)
      throw std::invalid_argument(std::string("ERROR: --trait1-sigsq must be non-negative"));
  }

  // initialize default value for --trait2-sigsq and check for out of range values
  if (simu_options.num_traits >= 2) {
    auto& target = simu_options.trait2_sigsq;
    if (!target.empty() && (target.size() != simu_options.num_components))
      throw std::invalid_argument(std::string("ERROR: Number of --trait2-sigsq values does not match the number of components"));
    if (target.empty())
      for (int i = 0; i < simu_options.num_components; i++) target.push_back(DEFAULT_SIGSQ);
    for (auto val: target) if (val < 0)
      throw std::invalid_argument(std::string("ERROR: --trait2-sigsq must be non-negative"));
  }

  // initialize default value for --rg and check for out of range values
  if (simu_options.num_traits >= 2) {
    auto& target = simu_options.rg;
    if (!target.empty() && (target.size() != simu_options.num_components))
      throw std::invalid_argument(std::string("ERROR: Number of --rg values does not match the number of components"));
    if (target.empty())
      for (int i = 0; i < simu_options.num_components; i++) target.push_back(DEFAULT_RG);
    for (auto val: target) if (val < -1.0f || val > 1.0f)
      throw std::invalid_argument(std::string("ERROR: --rg values must be between -1.0 and 1.0"));
  }

  // Read bfiles (one or 22)
  for (auto bfile: simu_options.bfiles) {
    log << "Opening " << bfile << "... ";
    auto pio_file = std::make_shared<PioFile>(bfile);
    pio_files->push_back(pio_file);
    log << "done, "
        << pio_file->num_samples() << " samples, "
        << pio_file->num_loci() << " variants.\n";
    if (simu_options.num_samples < 0) simu_options.num_samples = pio_file->num_samples();
    if (simu_options.num_samples != pio_file->num_samples())
      throw std::runtime_error("ERROR: Inconsistent number of subjects in --bfile-chr input");
    simu_options.num_variants += pio_file->num_loci();
  }
  simu_options.num_variants = pio_files->num_variants();
  if (simu_options.bfiles.size() > 1)
    log << "Found " << simu_options.num_variants << " variants in total.\n";

  if (!simu_options.causal_variants.empty() || !simu_options.causal_regions.empty()) {
    // the PioFiles::snp_to_id_ map will be used only if we need to look up variants by their RS name,
    // e.i. if user specified --causal-variants or --causal-regions flags.
    log << "Validate that all variants in --bfile / --bfile-chr have unique names (required for --causal-variants and --causal-regions options)...\n";
    pio_files->find_variant_index_map();
  }

  // Initialize or validate ncas/ncon parameters
  if (simu_options.cc) {
    if (simu_options.ncas.empty() || (simu_options.ncon.empty())) {
      for (int i = 0; i < simu_options.num_traits; i++) {
        simu_options.ncas.push_back(static_cast<int>(simu_options.k[i] * simu_options.num_samples));
        simu_options.ncon.push_back(simu_options.num_samples - simu_options.ncas[i]);
      }
    }
    for (int i = 0; i < simu_options.num_traits; i++) {
      int max_cas = static_cast<int>(simu_options.k[i] * simu_options.num_samples);
      int max_con = simu_options.num_samples - max_cas;
      if ((simu_options.ncas[i] <= 0) || (simu_options.ncas[i] > max_cas))
        throw std::invalid_argument(std::string("ERROR: --ncas out of range, must be 0 to k*N, where k is prevalence"));
      if ((simu_options.ncon[i] <= 0) || (simu_options.ncon[i] > max_con))
        throw std::invalid_argument(std::string("ERROR: --ncon out of range, must be 0 to (1-k)*N, where k is prevalence"));
    }
  }

  if (simu_options.causal_regions.empty()) {
    if (!simu_options.causal_pi.empty()) {
      for (int i = 0; i < simu_options.num_components; i++)
        simu_options.causal_n.push_back(static_cast<int>(simu_options.causal_pi[i] * simu_options.num_variants));
      simu_options.causal_pi.clear();
    }

    if (!simu_options.causal_n.empty()) {
      int num_causal_variants = 0;
      for (int i = 0; i < simu_options.num_components; i++) num_causal_variants += simu_options.causal_n[i];
      if (num_causal_variants > simu_options.num_variants)
        throw std::invalid_argument("ERROR: sum of --causal-pi must be less than 1.0, or sum of --causal-n must be less than number of variants present in the input file.");
    }
  } else {
    std::vector<int> region_per_variant;  
    for (int i = 0; i < simu_options.num_variants; i++) region_per_variant.push_back(-1);

    // Validate that --causal-regions have no overlap
    for (int component_index = 0; component_index < simu_options.num_components; component_index++) {
      std::vector<int> variant_indices;
      read_causal_variants_file(simu_options.causal_regions[component_index], &log, pio_files, &variant_indices, nullptr, nullptr);
      for (int variant_index: variant_indices) {
        if (region_per_variant[variant_index] != -1) {
          std::stringstream ss; ss << "ERROR: Variant '" << pio_files->get_locus(variant_index)->name << "' is duplicated in --causal-regions input; --causal-regions must be non-overlapping.";
          throw std::runtime_error(ss.str());
        }
        region_per_variant[variant_index] = component_index;
      }

      if (!simu_options.causal_pi.empty()) {
        simu_options.causal_n.push_back(static_cast<int>(simu_options.causal_pi[component_index] * variant_indices.size()));
        if (simu_options.causal_n.back() <= 0) {
          std::stringstream ss; ss << "--causal-pi " << simu_options.causal_pi[component_index] 
                                   << " is too low because --causal-region " << simu_options.causal_regions[component_index]
                                   << " contains only " <<  variant_indices.size() << " variants.";
          throw std::runtime_error(ss.str()); 
        }
      }

      if (simu_options.causal_n[component_index] > variant_indices.size()) {
        std::stringstream ss; ss << "--causal-n " << simu_options.causal_n[component_index] 
                                  << " is too high because --causal-region " << simu_options.causal_regions[component_index]
                                  << " contains only " <<  variant_indices.size() << " variants.";
        throw std::runtime_error(ss.str()); 
      }
    }
    simu_options.causal_pi.clear();
  }

  // if (vm.count("seed") == 0)
  //  log << "--seed option was set to " << simu_options.seed << "\n";
}

void describe_simu_options(SimuOptions& s, Logger& log) {
  log << "Options in effect (after applying default setting to non-specified parameters):\n";
  if (!s.bfile.empty()) log << "\t--bfile " << s.bfile << " \\\n";
  if (!s.bfile_chr.empty()) log << "\t--bfile-chr " << s.bfile_chr << " \\\n";
  if (s.qt) log << "\t--qt \\\n";
  if (s.cc) log << "\t--cc \\\n";
  log << "\t--num-traits " << s.num_traits << " \\\n";
  if (!s.k.empty()) log << "\t--k " << s.k << " \\\n";
  if (!s.ncas.empty()) log << "\t--ncas " << s.ncas << " \\\n";
  if (!s.ncon.empty()) log << "\t--ncon " << s.ncon << " \\\n";
  if (!s.hsq.empty()) log << "\t--hsq " << s.hsq << " \\\n";
  log << "\t--num-components " << s.num_components << " \\\n";
  if (!s.causal_variants.empty()) log << "\t--causal-variants " << s.causal_variants << " \\\n";
  if (!s.causal_n.empty()) log << "\t--causal-n " << s.causal_n << " \\\n";
  if (!s.causal_pi.empty()) log << "\t--causal-pi " << s.causal_pi << " \\\n";
  if (!s.causal_regions.empty()) log << "\t--causal-regions " << s.causal_regions << " \\\n";
  if (!s.trait1_sigsq.empty()) log << "\t--trait1-sigsq " << s.trait1_sigsq << " \\\n";
  if (!s.trait2_sigsq.empty()) log << "\t--trait2-sigsq " << s.trait2_sigsq << " \\\n";
  if (!s.rg.empty()) log << "\t--rg " << s.rg << " \\\n";
  if (!s.out.empty()) log << "\t--out " << s.out << " \\\n";
  log << "\t--seed " << s.seed << " \\\n";
  if (s.verbose) log << "\t--verbose \\\n";
  if (s.gcta_sigma) log << "\t--gcta-sigma \\\n";
  if (s.norm_effect) log << "\t--norm-effect \\\n";
  log << "\n";
}

void find_causals(const SimuOptions& simu_options, boost::mt19937& rng, Logger& log, PioFiles* pio_files,
                  std::vector<int>* component_per_variant,
                  std::vector<double>* effect1_per_variant,
                  std::vector<double>* effect2_per_variant) {
  // generate vector of length simu_options.num_varinats
  // elements in the vector can be
  //    -1 means "not causal in either component", or 
  //    0..(num_components-1) gives index of causal component

  component_per_variant->clear();  
  for (int i = 0; i < simu_options.num_variants; i++) component_per_variant->push_back(-1);

  if (!simu_options.causal_variants.empty()) {
    // user explicitly specified causal variants
    bool effect_sizes_from_file = false;   // whether variants had effect size provided, or not

    for (int component_index = 0; component_index < simu_options.num_components; ++component_index) {
      std::vector<int> variant_indices;
      std::vector<double> causal_effect1, causal_effect2;
      read_causal_variants_file(simu_options.causal_variants[component_index], &log, pio_files, &variant_indices, &causal_effect1, &causal_effect2);

      if (component_index == 0) {
        if (!causal_effect1.empty()) for (int i = 0; i < simu_options.num_variants; i++) effect1_per_variant->push_back(0.0);
        if (!causal_effect2.empty()) for (int i = 0; i < simu_options.num_variants; i++) effect2_per_variant->push_back(0.0);
        effect_sizes_from_file = !causal_effect1.empty();

        if (effect_sizes_from_file && causal_effect2.empty() && (simu_options.num_traits == 2)) {
          std::stringstream ss; ss << "ERROR: issue with --causal-variants files; effect sizes are not provided for the second trait, file '" << simu_options.causal_variants[0] << "'";
          throw std::runtime_error(ss.str());
        }

        if (effect_sizes_from_file && !causal_effect2.empty() && (simu_options.num_traits == 1)) {
          log << "WARNING: minor issue with --causal-variants files; effect sizes provided for the second trait will be ignored, '" << simu_options.causal_variants[0] << "'\n";
        }
      }

      if (effect_sizes_from_file && causal_effect1.empty()) {
        std::stringstream ss; ss << "ERROR: issue with --causal-variants files; effect sizes are present in '" << simu_options.causal_variants[0] << "' but not in '" << simu_options.causal_variants[component_index] << "'";
        throw std::runtime_error(ss.str());
      }

      for (int i = 0; i < variant_indices.size(); i++) {
        int variant_index = variant_indices[i];
        if (component_per_variant->at(variant_index) != -1) {
          std::stringstream ss; ss << "ERROR: Variant '" << pio_files->get_locus(variant_index)->name << "' is duplicated in --causal-variants input";
          throw std::runtime_error(ss.str());
        }

        component_per_variant->at(variant_index) = component_index;
        if (effect_sizes_from_file && !causal_effect1.empty()) effect1_per_variant->at(variant_index) = causal_effect1[i];
        if (effect_sizes_from_file && !causal_effect2.empty()) effect2_per_variant->at(variant_index) = causal_effect2[i];
      }
    }
  } else {
    // user requested to randomly distribute causal variants
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > random_integer(rng, boost::uniform_int<>());
      
    if (!simu_options.causal_regions.empty()) {
      // draw causal variants from user-defined regions
      for (int component_index = 0; component_index < simu_options.num_components; component_index++) {
        std::vector<int> variant_indices;
        read_causal_variants_file(simu_options.causal_regions[component_index], nullptr, pio_files, &variant_indices, nullptr, nullptr);
        std::random_shuffle(variant_indices.begin(), variant_indices.end(), random_integer);

        for (int causal_index = 0; causal_index < simu_options.causal_n[component_index]; causal_index++) {
          component_per_variant->at(variant_indices[causal_index]) = component_index;
        }
      }
    } else {
      // draw causal variants randomly across all variants in --bfile / --bfile-chr
      std::vector<int> variant_indices_permutation;
      for (int i = 0; i < simu_options.num_variants; i++) variant_indices_permutation.push_back(i);
      std::random_shuffle(variant_indices_permutation.begin(), variant_indices_permutation.end(), random_integer);

      int permuted_index = 0;
      for (int component_index = 0; component_index < simu_options.num_components; component_index++) {
        for (int causal_index = 0; causal_index < simu_options.causal_n[component_index]; causal_index++, permuted_index++) {
          component_per_variant->at(variant_indices_permutation[permuted_index]) = component_index;
        }
      }
    }
  }
}

void find_freq(const SimuOptions& simu_options, 
               const std::vector<int>& component_per_variant,
               PioFiles* pio_files, 
               std::vector<double>* freq_vec) {
  freq_vec->clear();
  std::vector<snp_t> buffer;
  for (int variant_index = 0; variant_index < simu_options.num_variants; variant_index++) {
    freq_vec->push_back(std::numeric_limits<double>::quiet_NaN());

    if (component_per_variant[variant_index] == -1) {
      // skip variant because we are not interested
      pio_status_t retval = pio_files->skip_row();
      if (retval == PIO_END) {
        std::stringstream ss; ss << "ERROR: expect to find " << simu_options.num_variants << " variants across input files; found only " << variant_index << ". Are input Plink files corrupted?";
        throw std::runtime_error(ss.str());
      }

      continue;
    }

    pio_status_t retval = pio_files->next_row(&buffer);
    if (retval == PIO_END) {
      std::stringstream ss; ss << "ERROR: expect to find " << simu_options.num_variants << " variants across input files; found only " << variant_index << ". Are input Plink files corrupted?";
      throw std::runtime_error(ss.str());
    }

    int genocount[4] = {0, 0, 0, 0};   // 0 = AA, 1 = Aa, 2 = aa, 3 = missing
    for (int sample_index = 0; sample_index < simu_options.num_samples; sample_index++) {
      genocount[(int)buffer[sample_index]]++;
    }

    int genocount_non_missing = genocount[2] + genocount[1] + genocount[0];
    freq_vec->at(variant_index) = (genocount_non_missing == 0) ? 0.0 : static_cast<double>(2 * genocount[0] + 1 * genocount[1]) / static_cast<double>(2 * genocount_non_missing);
  }

  // sanity-check -- now we should have reached end of the stream
  if (pio_files->skip_row() != PIO_END) {
    std::stringstream ss; ss << "ERROR: expect to find " << simu_options.num_variants << " variants across input files; more variants found. Are input Plink files corrupted?";
    throw std::runtime_error(ss.str());
  }

  pio_files->reset_row();
}

void find_effect_sizes(const SimuOptions& simu_options, boost::mt19937& rng,
                       const std::vector<int>& component_per_variant,
                       std::vector<double>* effect1_per_variant) {
  // generate vector of effect sizes (for trait1) for each variant
  effect1_per_variant->clear();
  for (int i = 0; i < simu_options.num_variants; i++) effect1_per_variant->push_back(0.0);

  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_normal(rng, boost::normal_distribution<>());

  for (int i = 0; i < simu_options.num_variants; i++) {
    if (component_per_variant[i] == -1) continue;
    double sigma = sqrt(simu_options.trait1_sigsq[component_per_variant[i]]);
    effect1_per_variant->at(i) = sigma * random_normal();
  }
}

void find_effect_sizes_bivariate(const SimuOptions& simu_options, boost::mt19937& rng,
                                 const std::vector<int>& component_per_variant,
                                 std::vector<double>* effect1_per_variant,
                                 std::vector<double>* effect2_per_variant) {
  // generate vector of effect sizes (for trait1 and trait2) for each variant
  effect1_per_variant->clear(); effect2_per_variant->clear();
  for (int i = 0; i < simu_options.num_variants; i++) {
    effect1_per_variant->push_back(0.0);
    effect2_per_variant->push_back(0.0);
  }

  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_normal(rng, boost::normal_distribution<>());

  // https://www2.stat.duke.edu/courses/Spring12/sta104.1/Lectures/Lec22.pdf
  for (int i = 0; i < simu_options.num_variants; i++) {
    if (component_per_variant[i] == -1) continue;
    double sigma1 = sqrt(simu_options.trait1_sigsq[component_per_variant[i]]);
    double sigma2 = sqrt(simu_options.trait2_sigsq[component_per_variant[i]]);
    double rg = simu_options.rg[component_per_variant[i]];
    double x1 = random_normal(); double x2 = random_normal();
    effect1_per_variant->at(i) = sigma1 * x1;
    effect2_per_variant->at(i) = sigma2 * (rg * x1 + sqrt(1-rg*rg) * x2);
  }
}

void find_pheno(const SimuOptions& simu_options, 
                const std::vector<int>& component_per_variant,
                const std::vector<double>& freq_vec,
                const std::vector<double>& effect1_per_variant,
                const std::vector<double>& effect2_per_variant,
                PioFiles* pio_files, 
                std::vector<double>* pheno1_per_sample,
                std::vector<double>* pheno2_per_sample) {
  const bool bivariate = (simu_options.num_traits == 2);
  std::vector<snp_t> buffer;

  pheno1_per_sample->clear(); pheno2_per_sample->clear();
  for (int sample_index = 0; sample_index < simu_options.num_samples; sample_index++) {
    pheno1_per_sample->push_back(0);
    if (bivariate) pheno2_per_sample->push_back(0);
  }

  for (int variant_index = 0; variant_index < simu_options.num_variants; variant_index++) {
    if (component_per_variant[variant_index] == -1) {
      // skip variant because we are not interested
      pio_status_t retval = pio_files->skip_row();
      if (retval == PIO_END) {
        std::stringstream ss; ss << "ERROR: expect to find " << simu_options.num_variants << " variants across input files; found only " << variant_index << ". Are input Plink files corrupted?";
        throw std::runtime_error(ss.str());
      }

      continue;
    }

    pio_status_t retval = pio_files->next_row(&buffer);
    if (retval == PIO_END) {
      std::stringstream ss; ss << "ERROR: expect to find " << simu_options.num_variants << " variants across input files; found only " << variant_index << ". Are input Plink files corrupted?";
      throw std::runtime_error(ss.str());
    }

    double average_A1_dosage = 2 * freq_vec[variant_index]; 
    for (int sample_index = 0; sample_index < simu_options.num_samples; sample_index++) {
      // missing value have 0 contribution to phenotype.
      // This is because best we can do is to replace missing value with "2x(allele frequency)" (e.i. average_A1_dosage),
      // which cancels out because we center by "2x(allele frequency)".
      if (buffer[sample_index] == 3) continue;
      double num_of_A1_copies = static_cast<double>(2 - buffer[sample_index]);  // 0 = AA, 1 = Aa, 2 = aa

      pheno1_per_sample->at(sample_index) += (num_of_A1_copies - average_A1_dosage) * effect1_per_variant[variant_index];
      if (bivariate) pheno2_per_sample->at(sample_index) += (num_of_A1_copies - average_A1_dosage) * effect2_per_variant[variant_index];
    }
  }

  // sanity-check -- now we should have reached end of the stream
  if (pio_files->skip_row() != PIO_END) {
    std::stringstream ss; ss << "ERROR: expect to find " << simu_options.num_variants << " variants across input files; more variants found. Are input Plink files corrupted?";
    throw std::runtime_error(ss.str());
  }

  pio_files->reset_row();
}

void apply_heritability(float hsq, boost::mt19937& rng, std::vector<double>* pheno_per_sample, std::vector<double>* effect_per_variant) {
  double pheno_var = 0.0;
  double n = static_cast<double>(pheno_per_sample->size());
  for (int i = 0; i < pheno_per_sample->size(); i++) pheno_var += (pheno_per_sample->at(i) * pheno_per_sample->at(i));
  pheno_var = pheno_var / n;

  std::vector<double> noise;
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > random_normal(rng, boost::normal_distribution<>());
  for (int i = 0; i < pheno_per_sample->size(); i++) noise.push_back(random_normal());
  
  double gen_coef = 0, env_coef = 0;
  if ((fabs(hsq) <= FLT_EPSILON) || (fabs(pheno_var) <= FLT_EPSILON)) {
    // This is an interesting case --- why there are some true casual SNPs, and yet heritability is 0?
    gen_coef = 0; env_coef = 1;
  } else {
    gen_coef = sqrt(hsq / pheno_var);
    env_coef = sqrt(1.0 - hsq);  // under assumption of var(noise)==1.
  }

  for (int i = 0; i < pheno_per_sample->size(); i++) pheno_per_sample->at(i) = pheno_per_sample->at(i) * gen_coef + noise[i] * env_coef;
  for (int i = 0; i < effect_per_variant->size(); i++) effect_per_variant->at(i) *= gen_coef;
}

void apply_liability_threshold(float k, int ncas, int ncon, boost::mt19937& rng, std::vector<double>* pheno_per_sample) {
  std::vector<int> index(pheno_per_sample->size(), 0);
  for (int i = 0 ; i != index.size() ; i++) index[i] = i;
  std::sort(index.begin(), index.end(), [&](const int& a, const int& b) { return (pheno_per_sample->at(a) < pheno_per_sample->at(b)); });
  // now, pheno[index[0]] is the smallest element of pheno.

  int n = static_cast<int>(pheno_per_sample->size());
  int max_ncas = floor(k * n);
  int max_ncon = n - max_ncas;

  std::vector<int> cas_inds, con_inds;
  for (int i = 0; i < max_ncon; i++) con_inds.push_back(i);
  for (int i = max_ncon; i < n; i++) cas_inds.push_back(i);
  
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > random_integer(rng, boost::uniform_int<>());
  std::random_shuffle(con_inds.begin(), con_inds.end(), random_integer);
  std::random_shuffle(cas_inds.begin(), cas_inds.end(), random_integer);

  if (con_inds.size() < ncon) throw(std::runtime_error("ERROR: too many controls requested"));
  if (cas_inds.size() < ncas) throw(std::runtime_error("ERROR: too many cases requested"));

  for (int i = 0; i < n; i++) pheno_per_sample->at(i) = -9;
  for (int i = 0; i < ncon; i++) pheno_per_sample->at(index[con_inds[i]]) = 1;
  for (int i = 0; i < ncas; i++) pheno_per_sample->at(index[cas_inds[i]]) = 2;
}

void save_pheno_file(const SimuOptions& simu_options, PioFiles* pio_files,
                     const std::vector<double>& pheno1_per_sample, const std::vector<double>& pheno2_per_sample) {
  std::ofstream file(simu_options.out_pheno);
  const bool bivariate = (simu_options.num_traits==2);
  
  file << "FID\tIID\ttrait1" << (bivariate ? "\ttrait2" : "") << "\n";
  for (int i = 0; i < pheno1_per_sample.size(); i++) {
    pio_sample_t* sample = pio_files->get_sample(i);
    file << sample->fid << "\t" << sample->iid << "\t";
    file << pheno1_per_sample[i];
    if (bivariate) file << "\t" << pheno2_per_sample[i];
    file << "\n";
  }
}

void save_causals_file(const SimuOptions& simu_options,
                       const std::string& out_causal, PioFiles* pio_files,
                       const std::vector<double>& freq_vec,
                       const std::vector<int>& component_per_variant,
                       const std::vector<double>& effect_per_variant) {
  std::ofstream file(out_causal);
  file << "SNP\tCHR\tPOS\tA1\tA2\tFRQ";
  for (int i = 0; i < simu_options.num_components; i++) file << "\tBETA_c" << (i+1);
  file << "\n";
  for (int variant_index = 0; variant_index < simu_options.num_variants; variant_index++) {
    if (component_per_variant[variant_index] != -1) {
      pio_locus_t* locus = pio_files->get_locus(variant_index);
      if (locus == nullptr) {
        std::stringstream ss; ss << "ERROR: unable to read locus with index " << variant_index;
        throw(std::runtime_error(ss.str()));
      }

      double effect = effect_per_variant[variant_index];
      if (simu_options.norm_effect) {
        double freq = freq_vec[variant_index];
        double het = fabs(2 * freq * (1-freq));
        effect *= sqrt(het);
      }

      file << locus->name << "\t";
      file << static_cast<int>(locus->chromosome) << "\t";
      file << locus->bp_position << "\t";
      file << locus->allele1 << "\t";
      file << locus->allele2 << "\t";
      file << freq_vec[variant_index];
      for (int i = 0; i < simu_options.num_components; i++)
        file << "\t" << ((i == component_per_variant[variant_index]) ? effect : 0.0);
      file << "\n";
    }
  }
}



int
main(int argc, char *argv[])
{
  try {
    SimuOptions simu_options;
    po::options_description po_options("SIMU " VERSION " - library for simulation of GWAS summary statistics");
    po_options.add_options()
      ("help,h", "produce this help message")    
      ("bfile", po::value(&simu_options.bfile), "prefix for plink .bed/.bim/.fam file")
      ("bfile-chr", po::value(&simu_options.bfile_chr),
        "same as --bfile, but will automatically concatenate .bed/.bim/.fam files split "
        "across 22 chromosomes. If the filename prefix contains the symbol @, SIMU will "
        "replace the @ symbol with chromosome numbers. Otherwise, SIMU will append chromosome "
        "numbers to the end of the filename prefix.")
      ("qt", po::bool_switch(&simu_options.qt)->default_value(false), "simulate quantitative trait")
      ("cc", po::bool_switch(&simu_options.cc)->default_value(false), "simulate case/control trait")
      ("num-traits", po::value(&simu_options.num_traits)->default_value(1), "number of traits (either 1 or 2 traits are supported)")
      ("k", po::value< std::vector<float> >(&simu_options.k)->multitoken(), "prevalence for case/control traits, by default 0.1; one value per trait")
      ("ncas", po::value< std::vector<int> >(&simu_options.ncas)->multitoken(), "number of cases, by default N*k, where N is sample size; one value per trait")
      ("ncon", po::value< std::vector<int> >(&simu_options.ncon)->multitoken(), "number of controls, by default N*(1-k), where N is sample size; one value per trait")
      ("hsq", po::value< std::vector<float> >(&simu_options.hsq)->multitoken(), "heritability, by default 0.7; one value per trait")
      ("num-components", po::value(&simu_options.num_components)->default_value(1), "number of components in the mixture")      
      ("causal-pi", po::value< std::vector<float> >(&simu_options.causal_pi)->multitoken(), "proportion of causal variants; by default 0.001; one value per mixture component")
      ("causal-n", po::value< std::vector<int> >(&simu_options.causal_n)->multitoken(), "number of causal variants (alternative to --causal-pi option); one value per mixture component")
      ("causal-variants", po::value< std::vector<std::string> >(&simu_options.causal_variants)->multitoken(),
       "file with a list of causal variants and, optionally, their effect sizes; "
       "one file per mixture component. "
       "This is an alternative option to --causal-pi and --causal-n. " 
       "See README.md file for detailed description of file formats.")
      ("causal-regions", po::value< std::vector<std::string> >(&simu_options.causal_regions)->multitoken(),
       "file with a list of non-overlapping regions to distribute causal variants. "
       "Regions must be defined as a list of variant names (e.g. RS numbers, one value per line). "
       "Supplements --causal-pi or --causal-n options; can not be used together with --causal-variants option. "
       "One file per mixture component")
      ("trait1-sigsq", po::value< std::vector<float> >(&simu_options.trait1_sigsq)->multitoken(), "variance of effect sizes for trait1 per causal marker; by default 1.0; one value per mixture component")
      ("trait2-sigsq", po::value< std::vector<float> >(&simu_options.trait2_sigsq)->multitoken(), "variance of effect sizes for trait2 per causal marker; by default 1.0; one value per mixture component")
      ("gcta-sigma", po::bool_switch(&simu_options.gcta_sigma)->default_value(false), "draw effect sizes with variance inversely proportional to sqrt(2*p(1-p)), where p is allele frequency.")
      ("norm-effect", po::bool_switch(&simu_options.norm_effect)->default_value(false),
       "report effect sizes w.r.t. normalized genotypes (e.i. additively coded 0,1,2 genotypes devided by sqrt(2*p(1-p)), where p is allele frequency). "
       "Default behavior without --norm-effect is to report effect size w.r.t. additively coded 0,1,2 genotypes.")
      ("rg", po::value< std::vector<float> >(&simu_options.rg)->multitoken(), "coefficient of genetic correlation; by default 0.0; one value per mixture component")
      ("seed", po::value(&simu_options.seed), "seed for random numbers generator (default is time-dependent seed)")
      ("out", po::value(&simu_options.out)->default_value("simu"),
       "prefix of the output files; will generate .pheno file containing synthesized phenotypes; "
      "and .*.causals files (one file per trait) containing lists of causal variants and their effect sizes for each component in the mixture. "
      "See README.md file for detailed description of file formats.")
      // ("verbose", po::bool_switch(&simu_options.verbose)->default_value(false), "enable verbose logging")
    ;

    // expand documentation - format of output files
    // for case-control traits, in pheno file, 1=unaffected (control), 2=affected (case); -9 for missing
    // feature - read causal effect sizes from file (instead of simulating them)
    // test and implement how to pass negative rg (since it starts with minus)

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .extra_style_parser(&ignore_numbers)
      .options(po_options)
      .run(), vm);
    notify(vm);

    bool show_help = (vm.count("help") > 0);
    if (show_help) {
      std::cerr << po_options;
      std::cerr << "\nExamples:\n";
      std::cerr << std::endl;
      std::cerr << "* Simulate a quantitative trait with 1% causal markers, with the heritability of 0.5:\n";
      std::cerr << "  simu --bfile test --qt --causal-pi 0.01 --hsq 0.5\n";
      std::cerr << std::endl;
      std::cerr << "* Simulate a quantitative trait as above, using simulation model from GCTA GWAS Simulation:\n";
      std::cerr << "  simu --bfile test --qt --causal-pi 0.01 --hsq 0.5 --gcta-sigma --norm-effect \n";
      std::cerr << std::endl;
      std::cerr << "* Simulate two quantitative traits with genetic correlation of 0.8 and the heritabilities of 0.2 and 0.6:\n";
      std::cerr << "  simu --bfile test --qt --causal-pi 0.01 --num-traits 2 --hsq 0.2 0.6 --rg 0.8 \n";
      std::cerr << std::endl;
      std::cerr << "* Simulate 500 cases and 500 controls with the heritability of liability of 0.5 and disease prevalence of 0.1:\n";
      std::cerr << "  simu --bfile test --cc --k 0.1 --ncas 500 --ncon 500 --causal-pi 0.01 --hsq 0.5 \n";
      std::cerr << std::endl;
      std::cerr << "* Simulate two quantitative traits with three genetic components:\n";
      std::cerr << "      (1) causal variants specific to the first trait, \n";
      std::cerr << "      (2) causal variants specific to the second trait, \n";
      std::cerr << "      (3) shared causal variants with correlation of 0.8: \n";
      std::cerr << "  simu --bfile test --qt --num-traits 2 --hsq 0.2 0.6  --num-components 3 \\ \n";
      std::cerr << "       --causal-pi 0.01 0.01 0.01 --trait1-sigsq 1 0 1 --trait2-sigsq 0 1 1 --rg 0 0 0.8 \n";
      std::cerr << std::endl;
      std::cerr << "* Simulate a quantitative trait with specific location of causal markers:\n";
      std::cerr << "  simu --bfile test --qt --causal-variants causal.snplist --hsq 0.5\n";
      std::cerr << std::endl;
      std::cerr << "* Simulate a quantitative trait where n=10 causal markers are randomly distributed within specific region:\n";
      std::cerr << "  simu --bfile test --qt --causal-n 10 --causal-regions snplist --hsq 0.5\n";
      std::cerr << std::endl;
      exit(EXIT_FAILURE);
    }

    Logger log(simu_options.out + ".log");
    log_header(argc, argv, log);

    try {
      auto analysis_started = boost::posix_time::second_clock::local_time();
      log << "Analysis started: " << analysis_started << "\n";

      // fix_and_validate needs to read input bfiles as it needs to know num_samples
      PioFiles pio_files;
      fix_and_validate(simu_options, vm, log, &pio_files);
      describe_simu_options(simu_options, log);

      // Generator of random numbers
      boost::mt19937 rng(simu_options.seed);
  
      // Find component for causal variants
      std::vector<int> component_per_variant;
      std::vector<double> effect1_per_variant, effect2_per_variant;
      find_causals(simu_options, rng, log, &pio_files, &component_per_variant, &effect1_per_variant, &effect2_per_variant);

      log << "Calculate allele frequencies (for causal markers only)... \n";
      std::vector<double> freq_vec;  // vector of allele frequency; NB! vector has -1 for non-causal variants, as we don't need them anywhere later;
      find_freq(simu_options, component_per_variant, &pio_files, &freq_vec);

      // Find effect sizes for causal variants (one per trait) --- unless user provided values in --causal-variants file
      if (effect1_per_variant.empty()) {
        log << "Generate effect sizes from normal distribution...\n";
        if (simu_options.num_traits==1) find_effect_sizes(simu_options, rng, component_per_variant, &effect1_per_variant);
        else find_effect_sizes_bivariate(simu_options, rng, component_per_variant, &effect1_per_variant, &effect2_per_variant);

        // Apply --gcta-sigma flag
        if (simu_options.gcta_sigma) {
          log << "Apply --gcta-sigma option to effect sizes...\n";
          for (int variant_index = 0; variant_index < simu_options.num_variants; variant_index++) {
            if (component_per_variant[variant_index] == -1) continue;
            double freq = freq_vec[variant_index];
            double het = fabs(2 * freq * (1-freq));
            if (het < FLT_EPSILON) {
              pio_locus_t* locus = pio_files.get_locus(variant_index);
              std::stringstream ss; ss << "Causal variant " << ((locus == nullptr) ? "rs???" : locus->name) << " has zero frequency. Can not apply --gcta-sigma.";
              throw std::runtime_error(ss.str());
            }

            effect1_per_variant[variant_index] /= sqrt(het);
            if (simu_options.num_traits==2) effect2_per_variant[variant_index] /= sqrt(het);
          }
        }
      } else {
        log << "Use effect sizes provided in --causal-variants file.\n";
        if (simu_options.num_components > 1)
          log << "WARNING: when --causal-variants contain effect sizes there is no real need to use more than one component (--num-components), isn't it so?\n";
        if (simu_options.gcta_sigma || (vm.count("trait1-sigsq") > 0) || (vm.count("trait2-sigsq") > 0) || (vm.count("rg") > 0))
          log << "WARNING: Options --gcta-sigma, --trait1-sigsq, --trait2-sigsq, --rg are irrelevant and will be ignored.\n";
        if (simu_options.norm_effect) {
          log << "NB! due to --norm-effect option, all effect sizes from --causal-variants will be treated as per-normalized-genotypes effect sizes\n";
          for (int variant_index = 0; variant_index < simu_options.num_variants; variant_index++) {
            if (component_per_variant[variant_index] == -1) continue;
            double freq = freq_vec[variant_index];
            double het = fabs(2 * freq * (1-freq));
            if (het < FLT_EPSILON) {
              pio_locus_t* locus = pio_files.get_locus(variant_index);
              std::stringstream ss; ss << "Causal variant " << ((locus == nullptr) ? "rs???" : locus->name) << " has zero frequency. Can not apply --norm-effect.";
              throw std::runtime_error(ss.str());
            }

            effect1_per_variant[variant_index] /= sqrt(het);
            if (simu_options.num_traits==2) effect2_per_variant[variant_index] /= sqrt(het);
          }
        } else {
          log << "NB! without --norm-effect option all effect sizes from --causal-variants will be treated w.r.t. additive coded 0,1,2 genotypes\n";
        }
      }

      log << "Calculate phenotypes... \n";
      // Find phenotypes
      std::vector<double> pheno1_per_sample, pheno2_per_sample;
      find_pheno(simu_options, component_per_variant, freq_vec, effect1_per_variant, effect2_per_variant,
                 &pio_files, &pheno1_per_sample, &pheno2_per_sample);

      // Apply heritability
      apply_heritability(simu_options.hsq[0], rng, &pheno1_per_sample, &effect1_per_variant);
      if (simu_options.num_traits==2) apply_heritability(simu_options.hsq[1], rng, &pheno2_per_sample, &effect2_per_variant);

      // Apply liability threshold model
      if (simu_options.cc) {
        apply_liability_threshold(simu_options.k[0], simu_options.ncas[0], simu_options.ncon[0], rng, &pheno1_per_sample);
        if (simu_options.num_traits==2) apply_liability_threshold(simu_options.k[1], simu_options.ncas[1], simu_options.ncon[1], rng, &pheno2_per_sample);
      }

      // Save phenotypes to output files
      log << "Save phenotypes to " << simu_options.out_pheno << "...\n";
      save_pheno_file(simu_options, &pio_files, pheno1_per_sample, pheno2_per_sample);

      // Save causals to output files
      log << "Save causal variants and their effect sizes to "
          << simu_options.out << ((simu_options.num_traits==2) ? ".*" : ".1") << ".causals...\n";
      save_causals_file(simu_options, simu_options.out_causals[0], &pio_files, freq_vec, component_per_variant, effect1_per_variant);
      if (simu_options.num_traits==2) save_causals_file(simu_options, simu_options.out_causals[1], &pio_files, freq_vec, component_per_variant, effect2_per_variant);

      if (simu_options.verbose) {  // debug printing
        log << "[VERBOSE] First 10 elements in resulting files:\n";
        log << "[VERBOSE]\t"; for (int i = 0; i < 10; i++) {if (i >= component_per_variant.size()) break; log << component_per_variant[i] << " "; }; log << "<- component\n";
        log << "[VERBOSE]\t"; for (int i = 0; i < 10; i++) {if (i >= freq_vec.size()) break; log << freq_vec[i] << " "; }; log << "<- freq(a1)\n";
        log << "[VERBOSE]\t"; for (int i = 0; i < 10; i++) {if (i >= effect1_per_variant.size()) break; log << effect1_per_variant[i] << " "; }; log << "<- effect1\n";
        log << "[VERBOSE]\t"; for (int i = 0; i < 10; i++) {if (i >= effect2_per_variant.size()) break; log << effect2_per_variant[i] << " "; }; log << "<- effect2\n";
        log << "[VERBOSE]\t"; for (int i = 0; i < 10; i++) {if (i >= pheno1_per_sample.size()) break; log << pheno1_per_sample[i] << " "; }; log << "<- pheno1\n";
        log << "[VERBOSE]\t"; for (int i = 0; i < 10; i++) {if (i >= pheno2_per_sample.size()) break; log << pheno2_per_sample[i] << " "; }; log << "<- pheno2\n";
      }

      auto analysis_finished = boost::posix_time::second_clock::local_time();
      log << "Analysis finished: " << analysis_finished << "\n";
      log << "Elapsed time: " << analysis_finished - analysis_started << "\n";
    } catch(std::exception& e) {
      log << e.what() << "\n";
      return EXIT_FAILURE;
    }
  } catch(std::exception& e) {
    std::cerr << "Exception  : " << e.what() << "\n";
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown error occurred.";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

