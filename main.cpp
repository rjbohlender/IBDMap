#include <iostream>
#include <boost/program_options.hpp>
#include <thread>
#include "src/parser.hpp"
#include "src/parameters.hpp"
#include "src/threadpool.hpp"
#include "src/statistic.hpp"
#include "src/reporter.hpp"
#include "src/geneticmap.hpp"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  // C++ IO only
  std::ios_base::sync_with_stdio(false);

  po::options_description desc;
  bool global = false;

  desc.add_options()
		  ("input,i",
		   po::value<std::string>()->required(),
		   "Unified IBD input file.")
		  ("pheno,p",
		   po::value<std::string>()->required(),
		   ".ped format file describing case-control status of all samples.")
		  ("average,a",
		   po::value<std::string>()->required(),
		   "The genome wide difference in average rate of ibd between pairs in the original data.")
		  ("gmap",
		   po::value<std::vector<std::string>>()->multitoken()->required(),
		   "Genetic map file locations for mean normalization.")
		  ("total_breakpoints",
		   po::value<double>()->required(),
		   "The total number of breakpoints across all chromosomes.")
		  ("threads,t",
		   po::value<unsigned long>()->default_value(std::thread::hardware_concurrency() / 2),
		   "Number of threads used in execution.")
		  ("permutations,n",
		   po::value<unsigned long>()->default_value(1000000),
		   "Number of permutations.")
		  ("successes,s",
		   po::value<unsigned long>()->default_value(30),
		   "Number of successful permutations for early termination.")
		  ("seed",
		   po::value<unsigned int>(),
		   "Specify the seed to be shared by all breakpoints for equal permutations.")
		  ("help,h", "Display this message.")
		  ("global,g", po::bool_switch(&global), "Run subtracting the genome wide average in each permutation.");

  po::variables_map vm;
  try {
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	  std::cerr << desc << std::endl;
	  return 1;
	}
	if (vm.count("global") > 0 && vm.count("gmap") == 0) {
	  std::cerr << desc << std::endl;
	  std::cerr << "ERROR: -g requires --gmap." << std::endl;
	  return 1;
	}
  } catch (po::error &pe) {
	std::cerr << pe.what() << std::endl;
	std::cerr << desc << std::endl;
	return 1;
  }
  std::cerr << "Running parser.\n";
  Parser<std::string> parser(vm["input"].as<std::string>(), vm["pheno"].as<std::string>(), vm["average"].as<std::string>());

  std::cerr << "Running parameters.\n";

  unsigned int seed;
  if (vm.count("seed") > 0) {
	seed = vm["seed"].as<unsigned int>();
  } else {
	seed = std::random_device{}();
  }

  Parameters parameters{
	  vm["permutations"].as<unsigned long>(),
	  vm["successes"].as<unsigned long>(),
	  vm["threads"].as<unsigned long>(),
	   seed,
	   vm["total_breakpoints"].as<double>()
  };

  if (parameters.nthreads < 3) {
	throw (std::runtime_error("ERROR: Need at least 3 threads."));
  }

  // Initialize reporter
  auto reporter = std::make_shared<Reporter>();
  if (vm["global"].as<bool>()) {
	unsigned long nperms = 0;
	ThreadPool<std::pair<double, double>> threadpool(parameters);
	std::vector<Statistic> stats;
	std::mt19937 gen(parameters.seed);
	for (unsigned long i = 0; i < parser.nbreakpoints; i++) {
	  Statistic stat(arma::sp_colvec(parser.data.col(i)),
					 parser.breakpoints[i],
					 parser.indexer,
					 parser,
					 reporter,
					 parameters);
	  stat.original = stat.calculate();
	  stats.emplace_back(stat);
	}

	auto shuffle = [&](std::vector<int> &p) {
	  for (int i = p.size() - 1; i > 0; i--) {
		std::uniform_int_distribution<> dis(0, i);
		int j = dis(gen);
		double tmp = p[j];
		p[j] = p[i];
		p[i] = tmp;
	  }
	};

	GeneticMap geneticMap(vm["gmap"].as<std::vector<std::string>>());

	arma::vec genetic_dist(parser.nbreakpoints, arma::fill::zeros);
	std::string prechr = "0";
	std::string chr = "0";
	int pos = 0;
	int prepos = 0;
	double predistance = 0;
	arma::uword current = 0;
	for (const auto &bp : parser.breakpoints) {
	  chr = bp.breakpoint.first;
	  if (!boost::starts_with(chr, "chr")) {
		std::stringstream ss;
		ss << "chr" << chr;
		chr = ss.str();
	  }
	  pos = std::stoi(bp.breakpoint.second);
	  std::cerr << "chr: " << chr << " pos: " << pos << std::endl;

	  auto variants = geneticMap.find_nearest(chr, pos);
	  std::cerr << "first: " << variants.first.first << " second: " << variants.second.first << std::endl;
	  double distance;
	  if (variants.first == variants.second
		  || variants.first.second == variants.second.second) { // Variant found or best we could do
		distance = variants.first.second;
	  } else { // Linear interpolation
		distance = variants.first.second
			+ static_cast<double>(pos - variants.first.first) / (variants.second.first - variants.first.first)
				* (variants.second.second - variants.first.second);
	  }
	  if (prechr == chr) {
		genetic_dist(current - 1) = distance - predistance;
		prechr = chr;
		prepos = pos;
		predistance = distance;
	  } else {
		prechr = chr;
		prepos = pos;
		predistance = distance;
	  }
	  current++;
	}

	std::vector<int> phenotypes = parser.phenotypes;
	// Something wrong with the initial counts. Not clear why it isn't working.
	std::cerr << "phenotypes: " << arma::conv_to<arma::vec>::from(phenotypes).t();
	arma::vec original(parser.nbreakpoints, arma::fill::zeros);
	arma::vec successes(parser.nbreakpoints, arma::fill::zeros);
	arma::vec permutations(parser.nbreakpoints, arma::fill::zeros);

	std::vector<std::future<std::pair<double, double>>> original_counts;

	for (unsigned long i = 0; i < parser.nbreakpoints; i++) {
	  std::packaged_task<std::pair<double, double>()> f(std::bind(&Statistic::calc_original_rates, &stats[i]));
	  original_counts.emplace_back(f.get_future());
	  threadpool.submit(std::move(f));
	}

	arma::vec cscs(parser.nbreakpoints, arma::fill::zeros); // case case rate
	arma::vec other(parser.nbreakpoints, arma::fill::zeros); // other rate -- case control or control control or both
	// Something wrong with the initial counts. Not clear why it isn't working.
	current = 0;
	for (auto &f: original_counts) {
	  auto res = f.get(); // Get counts for variant
	  cscs(current) = res.first;
	  other(current) = res.second;
	  current++;
	}

	double cscs_mean = arma::accu(cscs % genetic_dist) / arma::accu(genetic_dist);
	double other_mean = arma::accu(other % genetic_dist) / arma::accu(genetic_dist);

	std::cerr << "original -- cscs_mean: " << cscs_mean << " cscn_mean: " << other_mean << std::endl;

	original = (cscs - cscs_mean) / parser.indexer.case_case - (other - other_mean) / parser.indexer.case_cont;

	std::cerr << "cscs: " << cscs.t();
	std::cerr << "cscn: " << other.t();
	std::cerr << "original: " << original.t();

	while (nperms < vm["permutations"].as<unsigned long>()) {
	  shuffle(phenotypes);
	  std::vector<std::future<std::pair<double, double>>>
		  results; // results require promises to be alive at time of get
	  for (unsigned long i = 0; i < parser.nbreakpoints; i++) {
		std::packaged_task<std::pair<double, double>()> f(std::bind(&Statistic::calc_rates, &stats[i], phenotypes));
		results.emplace_back(f.get_future());
		threadpool.submit(std::move(f));
	  }

	  cscs = arma::vec(parser.nbreakpoints, arma::fill::zeros);
	  other = arma::vec(parser.nbreakpoints, arma::fill::zeros);
	  arma::vec diffs(parser.nbreakpoints, arma::fill::zeros);
	  current = 0;
	  for (auto &f: results) {
		auto res = f.get(); // Get counts for variant
		cscs(current) = res.first;
		other(current) = res.second;
		current++;
	  }

	  cscs_mean = arma::accu(cscs % genetic_dist) / arma::accu(genetic_dist);
	  other_mean = arma::accu(other % genetic_dist) / arma::accu(genetic_dist);

	  diffs = (cscs - cscs_mean) / parser.indexer.case_case - (other - other_mean) / parser.indexer.case_cont;
	  // Genome wide average difference in permutation is: cscs_d - other_d
	  for (arma::uword i = 0; i < diffs.n_elem; i++) {
		if (diffs(i) > original(i))
		  stats[i].successes++;
		stats[i].permutations++;
	  }
	  nperms++;
	}

	for (int i = 0; i < stats.size(); i++) {
	  std::stringstream iss;
	  iss << stats[i].bp.breakpoint.first << "\t" << stats[i].bp.breakpoint.second << "\tP-value:\t"
		  << (stats[i].successes + 1) / (stats[i].permutations + 1) << "\tSuccesses:\t" << stats[i].successes
		  << "\tPerms:\t" << stats[i].permutations
		  << "\tDelta:\t" << stats[i].original << "\tWeightedDelta:\t" << original(i) << "\tlength:\t"
		  << genetic_dist(i) << "\tcscs:\t" << stats[i].case_case_ibd << "\tcscn:\t"
		  << stats[i].case_cont_ibd << "\tcncn:\t"
		  << stats[i].cont_cont_ibd << std::endl;
	  reporter->submit(iss.str());
	}

  } else {
	// Initialize ThreadPool
	unsigned long submitted = 0;
	ThreadPool<void> threadpool(parameters);
	std::vector<Statistic> stats;
	stats.reserve(parser.nbreakpoints);
	for (unsigned long i = 0; i < parser.nbreakpoints; i++) {
	  Statistic stat(arma::sp_colvec(parser.data.col(i)),
					 parser.breakpoints[i],
					 parser.indexer,
					 parser,
					 reporter,
					 parameters);
	  // stat.run();
	  stats.emplace_back(std::move(stat));
	  std::packaged_task<void()> f(std::bind(&Statistic::run, &stats.back()));
	  threadpool.submit(std::move(f));
	  submitted++;
	}
	std::cerr << "nsubmitted: " << submitted << std::endl;
#if 1
	while (!std::all_of(stats.begin(), stats.end(), [](Statistic &s) { return s.done; })) {
	  std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
	}
#endif
  }

#if 0
  std::ofstream statout("stat.out.txt");
  std::ofstream cscsout("cscs.out.txt");
  std::ofstream cscnout("cscn.out.txt");
  for(const auto &stat : stats) {
	statout << stat.bp.breakpoint.first << "\t" << stat.bp.breakpoint.second;
	cscsout << stat.bp.breakpoint.first << "\t" << stat.bp.breakpoint.second;
	cscnout << stat.bp.breakpoint.first << "\t" << stat.bp.breakpoint.second;
	for(int i = 0; i < stat.permuted.size(); i++) {
	  statout << "\t" << stat.permuted[i];
	  cscsout << "\t" << stat.permuted_cscs[i];
	  cscnout << "\t" << stat.permuted_cscn[i];
	}
	statout << std::endl;
	cscsout << std::endl;
	cscnout << std::endl;
  }
#endif
}