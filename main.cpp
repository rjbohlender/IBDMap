#include <iostream>
#include <optional>
#include <thread>
#include <fmt-7.0.3/include/fmt/ostream.h>
#include "src/parser.hpp"
#include "src/parameters.hpp"
#include "src/reporter.hpp"
#include "src/geneticmap.hpp"
#include "CLI11.hpp"

int main(int argc, char *argv[]) {
  // C++ IO only
  std::ios_base::sync_with_stdio(false);

  arma::wall_clock timer;
  timer.tic();

  Parameters params;
  CLI::App app{"carvaIBD is an IBD mapping tool for large scale IBD datasets."};

  app.add_option("-i,--input", params.input, "Unified IBD input file.")->required();
  app.add_option("-p,--pheno",
				 params.pheno,
				 "Path to file containing sample phenotype pairs. 1 for affected, 0 for unaffected, NA for sample to be skipped. Header line is required.")
	 ->required();
  app.add_option("-g,--gmap",
				 params.gmap,
				 "Recombination map files. Assumed to be three columns, position, chromosome, cM. Header line is required.")
	 ->required();
  app.add_option("-c,--cov", params.cov, "Path to covariates. Expected format is ID Value1 Value2 ...");
  app.add_option("--info", params.info, "Path to supporting info file. Expected format is chr segID ...");
  app.add_option("-t,--threads", params.nthreads, "Number of threads used in execution.");
  app.add_option("-n,--permutations", params.nperms, "Number of permutations.");
  app.add_option("--lower_bound",
				 params.lower_bound,
				 "If set, establishes a lower bound on the number of pairs that must be present at breakpoint for it to be included.");
  app.add_option("-m,--min_dist",
				 params.min_dist,
				 "Sets the minimum genetic distance between sites. The parser will automatically skip breakpoints that are closer than the given distance. Default value of 0.0 cM.")
	 ->default_val(0.0);
  app.add_option("--cm",
				 params.cM,
				 "Sets the minimum segment length for segments to be included. The parser will automatically skip segments that are smaller than the given length.");
  app.add_option("--af",
				 params.AF,
				 "Sets the maximum allele frequency for segments to be included. The parser will automatically skip segments that are more frequent than the given cutoff.");
  app.add_option("-r,--rsquared",
				 params.rsquared,
				 "Sets the maximum correlation between sites. The parser will automatically skip breakpoints that are closer than the given distance. Default value of 1.0.");
  app.add_option("--range",
				 params.range,
				 "Range of positions, separated by a comma. e.g., 100000,250000 will only analyze breakpoints with positions in that range, inclusive of the endpoints.")
	 ->delimiter(',')->expected(2);
  app.add_option("--seed", params.seed, "Specify the seed to be shared by all breakpoints for equal permutations.");
  app.add_option("-o,--output", params.output_path, "Output to a specified file. Default output is stdout.");
  app.add_flag("--swap_pheno",
			   params.swap,
			   "Swap case-control status of individuals. Useful for examining problematic datasets.");
  app.add_flag("--contcont", params.contcont, "Use control/control pairs as part of statistic calculation.");
  app.add_flag("-v,--verbose", params.verbose, "Addtional diagnostic output.");
  app.add_flag("--enable_testing", params.enable_testing, "Enable testing functions for Statistic class.");

  CLI11_PARSE(app, argc, argv);

  fmt::print(std::cerr, "Input: {}\n", params.input);
  fmt::print(std::cerr, "Phenotypes: {}\n", params.pheno);
  fmt::print(std::cerr, "Gmap: {}\n", params.gmap);
  fmt::print(std::cerr, "Output: {}\n", params.output_path);
  if (params.cov) {
	fmt::print(std::cerr, "COV: {}\n", *params.cov);
  }
  if (params.info) {
	fmt::print(std::cerr, "INFO: {}\n", *params.info);
  }

  if (params.verbose) {
	fmt::print(std::cerr, "Reading genetic map.\n");
  }
  GeneticMap gmap(params.gmap);

  if (params.verbose) {
	fmt::print(std::cerr, "Constructing parameters.\n");
  }

  if (params.nthreads < 3) {
	throw (std::runtime_error("ERROR: Need at least 3 threads."));
  }
  // Initialize reporter
  auto reporter = std::make_shared<Reporter>(params.output_path);

  if (params.verbose) {
	fmt::print(std::cerr, "Running parser.\n");
  }
  Parser<std::string> parser(params.input,
							 params.pheno,
							 params.cov,
							 params,
							 reporter,
							 gmap);

  // Sort output
  reporter.reset();  // Ensure output is complete
  struct OutContainer {
	std::string chrom;
	int pos;
	std::vector<std::vector<std::string>> data;
  };
  if (!params.output_path.empty()) {
	std::ifstream reinput(params.output_path);

	std::string line;
	std::vector<OutContainer> sortable_output;
	while (std::getline(reinput, line)) {
	  RJBUtil::Splitter<std::string> splitter(line, " \t");
	  std::vector<std::string> stats;
	  if (splitter.size() > 0) {
		if (sortable_output.empty() || sortable_output.back().pos != std::stoi(splitter[1])) {
		  for (int i = 2; i < splitter.size(); i++) {
			stats.push_back(splitter[i]);
		  }
		  OutContainer oc{
			  splitter[0],
			  std::stoi(splitter[1]),
			  std::vector<std::vector<std::string>>()
		  };
		  oc.data.emplace_back(std::move(stats));
		  sortable_output.emplace_back(std::move(oc));
		} else {
		  // Continuing the output
		  for (int i = 2; i < splitter.size(); i++) {
			stats.push_back(splitter[i]);
		  }
		  sortable_output.back().data.emplace_back(std::move(stats));
		}
	  }
	}
	std::sort(sortable_output.begin(),
			  sortable_output.end(),
			  [](OutContainer &a, OutContainer &b) { return a.pos < b.pos; });
	reinput.close();
	std::ofstream reoutput(params.output_path);
	for (const auto &v: sortable_output) {
	  for (const auto &w : v.data) {
		fmt::print(reoutput, "{}\t{}\t{}\n", v.chrom, v.pos, fmt::join(w.begin(), w.end(), "\t"));
	  }
	}
	reoutput.close();
  }
  fmt::print(std::cerr, "Total runtime: {}\n", timer.toc());
}