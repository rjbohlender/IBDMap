#include "CLI11.hpp"
#include "src/geneticmap.hpp"
#include "src/parameters.hpp"
#include "src/parser.hpp"
#include "src/phenotypes.hpp"
#include "src/reporter.hpp"
#include <fmt/include/fmt/ostream.h>
#include <iostream>
#include <optional>
#include <sys/stat.h>

// 32bit xorshift PRNG
static inline uint32_t prng_u32(uint32_t p)
{
    uint32_t  state = p;
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

int main(int argc, char *argv[]) {
  // C++ IO only
  std::ios_base::sync_with_stdio(false);

  arma::wall_clock timer;
  timer.tic();

  Parameters params;
  std::vector<int> tmp_range_vec;
  std::vector<std::string> tmp_exclude_vec;
  std::vector<std::string> tmp_sample_vec;

  CLI::App app{"IBDMap is an IBD mapping tool for large scale IBD datasets."};

  app.add_option("-i,--input",
				 params.input,
				 "Unified IBD input file.")->required()->check(CLI::ExistingFile);
  app.add_option("-p,--pheno",
				 params.pheno,
				 "Path to file containing sample phenotype pairs. 1 for "
				 "affected, 0 for unaffected, NA for sample to be skipped. "
				 "Header line is required.")->required()->check(CLI::ExistingFile);
  app.add_option("-g,--gmap",
				 params.gmap,
				 "Recombination map files. Assumed to be three columns, "
				 "position, chromosome, cM. Header line is required. "
	 			 "Expected format is physical_pos\tchromosome\tgenetic_position")->required()->check(CLI::ExistingFile);
  app.add_option("--info",
				 params.info,
				 "Path to supporting info file. Expected format is chr segID ...")->check(CLI::ExistingFile);
  app.add_option("--cov",
                 params.cov,
                 "Path to supporting covariates file. Expected format is sampleID cov1 cov2 ...")->check(CLI::ExistingFile);
  app.add_option("-t,--threads",
				 params.nthreads,
				 "Number of threads used in execution. 2 threads are reserved for parsing and output.")->default_val(std::thread::hardware_concurrency() - 1);
  app.add_option("-n,--permutations",
				 params.nperms,
				 "Number of permutations for each breakpoint.")->default_val(100000);
  app.add_option("--lower_bound",
				 params.lower_bound,
				 "If set, establishes a lower bound on the number of pairs "
				 "that must be present at breakpoint for it to be included.");
  app.add_option("-m,--min_dist",
				 params.min_dist,
				 "Sets the minimum genetic distance between breakpoints. "
				 "The parser will automatically skip breakpoints that are "
				 "closer than the given distance. Default value of 0.05 cM.")->default_val(0.05);
  app.add_option("--cm",
				 params.cM,
				 "Sets the minimum segment length for segments to be included. "
				 "The parser will automatically skip segments that are smaller "
				 "than the given length. Default value of 3.0.")->default_val(3.0);
  app.add_option("--af",
				 params.AF,
				 "Sets the maximum allele frequency for segments to be "
				 "included. The parser will automatically skip segments that "
				 "are more frequent than the given cutoff. Default value of 0.01")->default_val(0.01);
  app.add_option("-r,--rsquared",
				 params.rsquared,
				 "Sets the maximum correlation between sites. The parser will "
				 "automatically skip breakpoints that are more correlated than the given "
				 "threshold.");
  app.add_option("--range",
				 tmp_range_vec,
				 "Range of positions, separated by a comma. e.g., 100000,250000 "
				 "will only analyze breakpoints with positions in that range, "
				 "inclusive of the endpoints.")
	 ->delimiter(',')->expected(2);
  app.add_option("--exclude",
               tmp_exclude_vec,
               "A series of ranges of the form 123-456,457-900 to exclude from the analyzed chromosome."
               "All breakpoints within the range (inclusive) will be dropped.")
        ->delimiter(',');
  app.add_option("--sample_list,-l",
				 tmp_sample_vec,
				 "Comma-separated list of samples to include in analysis. All "
				 "pairs not including one of these samples will be dropped.")
	 ->delimiter(',');
  app.add_option("--seed",
				 params.seed,
				 "Specify the seed to be shared by all breakpoints for "
				 "equal permutations.");
  app.add_option("--threshold",
                 params.threshold,
                 "The number of ibd segments that a breakpoint must have to be tested.")->default_val(0);
  app.add_option("-o,--output",
				 params.output_path,
				 "Output to a specified file. Default output is is ./output/{date-time}.results.zst.");
  app.add_option("--pheno_col",
                 params.pheno_col,
                 "The column to read from the phenotype file. Column 0 should be sample ID and column 1+ should be 0, 1, NA.")->default_val(1);
  app.add_option("--read_permutations",
                 params.read_permutations,
                 "Read the set of permutations from a pre-generated TSV file at the given path.")->check(CLI::ExistingFile);
  app.add_option("--epsilon",
                 params.epsilon,
                 "Set the value of epsilon for bin creation in covariate adjusted permutation. Larger values increase bin size and reduce number of bins, improving speed with loss of accuracy.")->default_val(0.0);
  app.add_option("--algorithm",
                 params.algorithm,
                 "Set the optimizer for fitting the GLM prior to covariate adjustment. Options are: irls, irls_svdnewton, irls_qr, irls_qr_R, gradient_descent.")->default_val("irls");
  app.add_flag("--swap_pheno",
			   params.swap,
			   "Swap case-control status of individuals. Useful for "
			   "examining problematic datasets.");
  app.add_flag("--contcont",
			   params.contcont,
			   "Use control/control pairs as part of statistic calculation.");
  app.add_flag("-v,--verbose",
			   params.verbose,
			   "Addtional diagnostic output.");
  app.add_flag("--enable_testing",
			   params.enable_testing,
			   "Enable testing functions for Statistic class.");
  app.add_flag("--dash",
			   params.dash,
			   "Expect input to be formatted by DASH.");
  app.add_flag("--print_debug",
               params.print_debug,
               "Print statistic component information to a text file.");
  app.add_flag("--old",
               params.oldformat,
               "Expect non-DASH input to include the npairs and nsegments columns.");
  app.add_flag("--compressed-memory",
               params.compressed_memory,
               "Use significantly less memory at the cost of 3x-5x longer runtime.");
  app.add_flag("--cscs_only",
               params.cscs_only,
               "Use only the case/case rate when outputting test statistics.");

  CLI11_PARSE(app, argc, argv);

  // Default output location
  if (params.output_path.empty()) {
    const char* out_dir = "./output";
    struct stat buffer;
    if(stat(out_dir, &buffer) != 0) {
      mkdir(out_dir, 0755);
	}
    char cur_time[100];
    std::time_t t = std::time(nullptr);
    std::strftime(cur_time, 99, "%F-%T", std::localtime(&t));
	std::stringstream default_output;
	default_output << out_dir << "/" << cur_time << ".results.zst";
	params.output_path = default_output.str();
  }

  // Have to handle this way because optional wrapped vector arguments don't seem to be supported.
  if (tmp_range_vec.size() == 2) {
	params.range = tmp_range_vec;
  }
  if (!tmp_exclude_vec.empty()) {
      std::vector<std::pair<int, int>> exclude_vec;
      for (auto &v : tmp_exclude_vec) {
          RJBUtil::Splitter<std::string> splitter(v, "-");
          exclude_vec.emplace_back(std::stoi(splitter[0]), std::stoi(splitter[1]));
      }
      params.exclude = exclude_vec;
  }
  if (!tmp_sample_vec.empty()) {
	params.sample_list = std::set<std::string>(tmp_sample_vec.begin(), tmp_sample_vec.end());
  }

  if (params.info && !params.dash) {
    fmt::print(std::cerr, "Info files are only supported with the --dash option.");
    std::exit(-1);
  }

  if (params.info && !params.dash) {
    fmt::print(std::cerr, "Info files are only supported with the --dash option.");
    std::exit(-1);
  }

  params.print(std::cerr);
  GeneticMap gmap(params.gmap);

  if (params.verbose) {
	fmt::print(std::cerr, "Constructing parameters.\n");
  }

  // Initialize reporter
  auto reporter = std::make_shared<Reporter>(params.output_path, params.print_debug);

  if (params.verbose) {
	fmt::print(std::cerr, "Running parser.\n");
  }
  std::array<uint32_t, 4> seed_vals = {0};
  seed_vals[0] = params.seed;
  for (int i = 1; i < 4; ++i) {
      seed_vals[i] = prng_u32(seed_vals[i - 1]);
  }
  std::seed_seq seed_source(seed_vals.begin(), seed_vals.end());

  if (params.compressed_memory) {
      Phenotypes<compressed_pheno_vector> pheno(params, seed_source);
      Parser<compressed_pheno_vector> parser(
              params,
              reporter,
              gmap,
              pheno);
  } else {
      Phenotypes<pheno_vector> pheno(params, seed_source);
      Parser<pheno_vector> parser(
              params,
              reporter,
              gmap,
              pheno);
  }

  // Sort output
  fmt::print(std::cerr, "Sorting output.\n");
  reporter->sort();
  reporter.reset();  // Ensure output is complete
  fmt::print(std::cerr, "Total runtime: {}\n", timer.toc());
}