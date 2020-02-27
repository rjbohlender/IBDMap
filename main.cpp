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
  bool swap = false;
  bool contcont = false;

  boost::optional<arma::uword> lower_bound;

  desc.add_options()
      ("input,i",
       po::value<std::string>()->required(),
       "Unified IBD input file.")
      ("pheno,p",
       po::value<std::string>()->required(),
       "Path to file containing sample phenotype pairs. 1 for affected, 0 for unaffected, NA for sample to be skipped. Header line is required.")
      ("gmap,g",
       po::value<std::vector<std::string>>()->required()->composing(),
       "Recombination map files. Assumed to be three columns, position, chromosome, cM. Header line is required.")
      ("cov,c",
       po::value<boost::optional<std::string>>()->default_value(boost::none, ""),
       ".ped format file describing case-control status of all samples.")
      ("threads,t",
       po::value<unsigned long>()->default_value(std::thread::hardware_concurrency() / 2),
       "Number of threads used in execution.")
      ("permutations,n",
       po::value<unsigned long>()->default_value(1000000),
       "Number of permutations.")
      ("lower_bound",
       po::value(&lower_bound),
       "If set, establishes a lower bound on the number of pairs that must be present at breakpoint for it to be included.")
      ("min_dist,m",
       po::value<double>()->default_value(0.0),
       "Sets the minimum genetic distance between sites. The parser will automatically skip breakpoints that are closer than the given distance. Default value of 0.0 cM.")
      ("seed",
       po::value<unsigned int>(),
       "Specify the seed to be shared by all breakpoints for equal permutations.")
      ("output,o",
       po::value<std::string>(),
       "Output to a specified file. Default output is stdout.")
      ("swap_pheno",
       po::bool_switch(&swap),
       "Swap case-control status of individuals. Useful for examining problematic datasets.")
      ("contcont",
       po::bool_switch(&contcont),
       "Use control/control pairs as part of statistic calculation.")
      ("help,h", "Display this message.");
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
  arma::wall_clock timer;
  timer.tic();

  if(swap) {
    std::cerr << "Swap phenotype activated.\n";
  }

  unsigned int seed;
  if (vm.count("seed") > 0) {
    seed = vm["seed"].as<unsigned int>();
  } else {
    seed = std::random_device{}();
  }

  std::cerr << "Reading genetic map.\n";
  GeneticMap gmap(vm["gmap"].as<std::vector<std::string>>());

  std::cerr << "Constructing parameters.\n";
  Parameters parameters{
      vm["permutations"].as<unsigned long>(),
      vm["threads"].as<unsigned long>(),
      vm.count("output") > 0 ? vm["output"].as<std::string>() : "",
      seed,
      lower_bound,
      swap,
      contcont,
      vm["min_dist"].as<double>()
  };

  if (parameters.nthreads < 3) {
    throw (std::runtime_error("ERROR: Need at least 3 threads."));
  }
  // Initialize reporter
  auto reporter = std::make_shared<Reporter>(parameters.output_path);

  std::cerr << "Running parser.\n";
  Parser<std::string> parser(vm["input"].as<std::string>(),
                             vm["pheno"].as<std::string>(),
                             vm["cov"].as<boost::optional<std::string>>(),
                             parameters,
                             reporter,
                             gmap);

  // Sort output
  reporter.reset();  // Ensure output is complete
  struct OutContainer {
    std::string chrom;
    int pos;
    std::vector<std::vector<std::string>> data;
  };
  if (!parameters.output_path.empty()) {
    std::ifstream reinput(parameters.output_path);

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
    std::ofstream reoutput(parameters.output_path);
    for (const auto &v: sortable_output) {
      for (const auto &w : v.data) {
        reoutput << v.chrom << "\t" << v.pos;
        for (const auto &s : w) {
          reoutput << "\t" << s;
        }
        reoutput << std::endl;
      }
    }
    reoutput.close();
  }
  std::cerr << "Total runtime: " << timer.toc() << std::endl;
}