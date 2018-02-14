#define VERSION "v0.9.0"

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
namespace po = boost::program_options;

#include <plinkio/plinkio.h>

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
  std::string keep;
  std::string out;

  // derived options
  std::vector<std::string> bfiles;
  int num_samples;
  std::string out_pheno;
  std::vector<std::string> out_causals;

  SimuOptions() {
    num_samples = -1;
  }
};

class Logger {
 public:
  Logger(std::string path) : log_file_(path) { }

  template <typename T>
  Logger& operator<< ( const T& rhs) {
    std::cout << rhs;
    log_file_ << rhs;
    return *this;
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

class PioFile {
 public:
  PioFile(std::string bfile) {
    if( pio_open( &plink_file_, bfile.c_str() ) != PIO_OK )
    {
      std::stringstream ss;
      ss << "ERROR: Could not open " << bfile;
      throw std::runtime_error(ss.str());
    }

    if( !pio_one_locus_per_row( &plink_file_ ) )
    {
      std::stringstream ss;
      ss << "ERROR: Unsupported format in " << bfile << ", snps must be rows, samples must be columns";
      throw std::runtime_error(ss.str());
    }

    row_size_ = pio_row_size( &plink_file_ );
    num_samples_ = pio_num_samples( &plink_file_ );
    num_loci_ = pio_num_loci( &plink_file_ );
  }

  int row_size() const { return row_size_; }
  int num_samples() const { return num_samples_; }
  int num_loci() const { return num_loci_; }

  ~PioFile() {
    pio_close( &plink_file_ );
  }

 private:
  pio_file_t plink_file_;
  int row_size_;
  int num_samples_;
  int num_loci_;
};

void fix_and_validate(SimuOptions& simu_options, po::variables_map& vm, Logger& log,
                      std::vector<std::shared_ptr<PioFile>>* pio_files) {
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
    simu_options.out_causals.push_back(simu_options.out + "." + boost::lexical_cast<std::string>(i) + ".causals");

  if ( boost::filesystem::exists( simu_options.out_pheno ) )
  {
    log << "WARNING: Target file " << simu_options.out_pheno << " already exists and will be overwritten\n"; 
  }

  if (simu_options.qt && (vm.count("k") > 0 || vm.count("ncas") > 0 || vm.count("ncon") > 0))
    log << "WARNING: Options --k, --ncas, --ncon are not relevant to --qt, and will be ignored\n";
  if ((simu_options.num_traits==1) && (vm.count("trait2-sigsq") > 0))
    log << "WARNING: Option --trait2-sigsq is not relevant for a single-trait simulations, and will be ignored\n";

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

  // initialize default value for k (prevalence for case/control trait) and check for out of range values
  if (simu_options.k.empty())
    for (int i = 0; i < simu_options.num_traits; i++) simu_options.k.push_back(0.1f);
  for (auto val: simu_options.k) if (val <= 0 || val >= 1)
    throw std::invalid_argument(std::string("ERROR: Prevalence --k must be between 0.0 and 1.0"));

  // inidialize default value for hsq (heritability) and check for out of range values
  if (simu_options.hsq.empty())
    for (int i = 0; i < simu_options.num_traits; i++) simu_options.hsq.push_back(0.7f);
  for (auto val: simu_options.hsq) if (val <= 0 || val >= 1)
    throw std::invalid_argument(std::string("ERROR: Heritability --hsq must be between 0.0 and 1.0"));

  // initialize default value for --causal-pi and check for out of range values
  if (!simu_options.causal_pi.empty() && (simu_options.causal_pi.size() != simu_options.num_components))
    throw std::invalid_argument(std::string("ERROR: Number of --causal-pi values does not match the number of components"));
  if (simu_options.causal_pi.empty())
    for (int i = 0; i < simu_options.num_components; i++) simu_options.causal_pi.push_back(0.001f);
  for (auto val: simu_options.causal_pi) if (val <= 0 || val >= 1)
    throw std::invalid_argument(std::string("ERROR: --causal-pi values must be between 0.0 and 1.0"));

  // initialize default value for --trait1-sigsq and check for out of range values
  {
    auto& target = simu_options.trait1_sigsq;
    if (!target.empty() && (target.size() != simu_options.num_components))
      throw std::invalid_argument(std::string("ERROR: Number of --trait1-sigsq values does not match the number of components"));
    if (target.empty())
      for (int i = 0; i < simu_options.num_components; i++) target.push_back(0.001f);
    for (auto val: target) if (val < 0)
      throw std::invalid_argument(std::string("ERROR: --trait1-sigsq must be non-negative"));
  }

  // initialize default value for --trait2-sigsq and check for out of range values
  if (simu_options.num_traits >= 2) {
    auto& target = simu_options.trait2_sigsq;
    if (!target.empty() && (target.size() != simu_options.num_components))
      throw std::invalid_argument(std::string("ERROR: Number of --trait2-sigsq values does not match the number of components"));
    if (target.empty())
      for (int i = 0; i < simu_options.num_components; i++) target.push_back(0.001f);
    for (auto val: target) if (val < 0)
      throw std::invalid_argument(std::string("ERROR: --trait2-sigsq must be non-negative"));
  }

  // initialize default value for --rg and check for out of range values
  {
    auto& target = simu_options.rg;
    if (!target.empty() && (target.size() != simu_options.num_components))
      throw std::invalid_argument(std::string("ERROR: Number of --rg values does not match the number of components"));
    if (target.empty())
      for (int i = 0; i < simu_options.num_components; i++) target.push_back(0.0f);
    for (auto val: target) if (val < -1.0f || val > 1.0f)
      throw std::invalid_argument(std::string("ERROR: --rg values must be between -1.0 and 1.0"));
  }

  // Read bfiles (one or 22)
  for (auto bfile: simu_options.bfiles) {
    log << "Opening " << bfile << "... ";
    pio_files->push_back(std::make_shared<PioFile>(bfile));
    log << "done, "
        << pio_files->back()->num_samples() << " samples, "
        << pio_files->back()->num_loci() << " variants.\n";
    if (simu_options.num_samples < 0) simu_options.num_samples = pio_files->back()->num_samples();
    if (simu_options.num_samples != pio_files->back()->num_samples())
      throw std::runtime_error("ERROR: Inconsistent number of subjects in --bfile-chr input");
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
}

int
main(int argc, char *argv[])
{
  try {
    SimuOptions simu_options;
    po::options_description po_options("SIMU " VERSION " - library for simulation of GWAS summary statistics");
    po_options.add_options()
      ("help,h", "produce help message")    
      ("bfile", po::value(&simu_options.bfile), "Prefix for Plink .bed/.bim/.fam file")
      ("bfile-chr", po::value(&simu_options.bfile_chr),
        "Same as --bfile, but will automatically concatenate .bed/.bim/.fam files split "
        "across 22 chromosomes. If the filename prefix contains the symbol @, SIMU will "
        "replace the @ symbol with chromosome numbers. Otherwise, SIMU will append chromosome "
        "numbers to the end of the filename prefix.")
      ("qt", po::bool_switch(&simu_options.qt)->default_value(false), "simulate quantitative trait")
      ("cc", po::bool_switch(&simu_options.cc)->default_value(false), "simulate case/control trait")
      ("num-traits", po::value(&simu_options.num_traits)->default_value(1), "Number of traits (either 1 or 2 traits are supported)")
      ("k", po::value< std::vector<float> >(&simu_options.k)->multitoken(), "prevalence for case/control traits, by default 0.1; one value per trait")
      ("ncas", po::value< std::vector<int> >(&simu_options.ncas)->multitoken(), "number of cases, by default N*k; one value per trait")
      ("ncon", po::value< std::vector<int> >(&simu_options.ncon)->multitoken(), "number of controls, by default N*(1-k); one value per trait")
      ("hsq", po::value< std::vector<float> >(&simu_options.hsq)->multitoken(), "heritability, by default 0.7; one value per trait")
      ("num-components", po::value(&simu_options.num_components)->default_value(1), "Number of components in the mixture")      
      ("causal-pi", po::value< std::vector<float> >(&simu_options.causal_pi)->multitoken(), "proportion of causal variants; by default 0.001; one value per mixture component")
      ("trait1-sigsq", po::value< std::vector<float> >(&simu_options.trait1_sigsq)->multitoken(), "variance of effect sizes for trait1 per causal marker; by default 1.0; one value per mixture component")
      ("trait2-sigsq", po::value< std::vector<float> >(&simu_options.trait2_sigsq)->multitoken(), "variance of effect sizes for trait2 per causal marker; by default 1.0; one value per mixture component")
      ("rg", po::value< std::vector<float> >(&simu_options.rg)->multitoken(), "coefficient of genetic correlation; by default 0.0; one value per mixture component")
      ("out", po::value(&simu_options.out)->default_value("simu"), "Prefix of the output file; will generate .pheno file (phenotypes) and .1.causals file (one per trait, list MarkerName for all causal variants and their effect sizes.")
    ;

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(po_options).run(), vm);
    notify(vm);
   
    bool show_help = (vm.count("help") > 0);
    if (show_help) {
      std::cerr << po_options;
      exit(EXIT_FAILURE);
    }

    Logger log(simu_options.out + ".log");
    log_header(argc, argv, log);

    try {
      // fix_and_validate needs to read input bfiles as it needs to know num_samples
      std::vector<std::shared_ptr<PioFile>> pio_files;
      fix_and_validate(simu_options, vm, log, &pio_files);

      // TBD - implement the logic
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

