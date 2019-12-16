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
      ("threads,t",
       po::value<unsigned long>()->default_value(std::thread::hardware_concurrency() / 2),
       "Number of threads used in execution.")
      ("permutations,n",
       po::value<unsigned long>()->default_value(1000000),
       "Number of permutations.")
      ("seed",
       po::value<unsigned int>(),
       "Specify the seed to be shared by all breakpoints for equal permutations.")
      ("output,o",
       po::value<std::string>(),
       "Output to a specified file. Default output is stdout.")
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

  unsigned int seed;
  if (vm.count("seed") > 0) {
    seed = vm["seed"].as<unsigned int>();
  } else {
    seed = std::random_device{}();
  }

  std::cerr << "Running parameters.\n";
  Parameters parameters{
      vm["permutations"].as<unsigned long>(),
      vm["threads"].as<unsigned long>(),
      vm.count("output") > 0 ? vm["output"].as<std::string>() : "",
      seed
  };

  if (parameters.nthreads < 3) {
    throw (std::runtime_error("ERROR: Need at least 3 threads."));
  }
  // Initialize reporter
  auto reporter = std::make_shared<Reporter>(parameters.output_path);

  std::cerr << "Running parser.\n";
  Parser<std::string> parser(vm["input"].as<std::string>(), vm["pheno"].as<std::string>(), parameters, reporter);

  // Sort output
  reporter.reset();  // Ensure output is complete
  struct OutContainer {
    std::string chrom;
    int pos;
    std::vector<std::string> data;
  };
  if (!parameters.output_path.empty()) {
    std::ifstream reinput(parameters.output_path);

    std::string line;
    std::vector<OutContainer> sortable_output;
    while(std::getline(reinput, line)) {
      RJBUtil::Splitter<std::string> splitter(line, " \t");
      std::vector<std::string> stats;
      if (splitter.size() > 0) {
        for (int i = 2; i < splitter.size(); i++) {
          stats.push_back(splitter[2]);
        }
        OutContainer oc {
            splitter[0],
            std::stoi(splitter[1]),
            std::move(stats)
        };
        sortable_output.emplace_back(std::move(oc));
      }
    }
    std::sort(sortable_output.begin(), sortable_output.end(), [](OutContainer &a, OutContainer &b){ return a.pos < b.pos; });
    reinput.close();
    std::ofstream reoutput(parameters.output_path);
    for (const auto &v: sortable_output) {
      reoutput << v.chrom << "\t" << v.pos;
      for (const auto &s : v.data) {
        reoutput << "\t" << s;
      }
      reoutput << std::endl;
    }
    reoutput.close();
  }
  std::cerr << "Total runtime: " << timer.toc() << std::endl;
}