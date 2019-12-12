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
            ("successes,s",
             po::value<unsigned long>()->default_value(30),
             "Number of successful permutations for early termination.")
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
    std::cerr << "Running parser.\n";
    Parser<std::string> parser(vm["input"].as<std::string>(), vm["pheno"].as<std::string>());

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
            vm.count("output") > 0 ? vm["output"].as<std::string>() : "",
            seed
    };

    if (parameters.nthreads < 3) {
        throw (std::runtime_error("ERROR: Need at least 3 threads."));
    }

    // Initialize reporter
    auto reporter = std::make_shared<Reporter>(parameters.output_path);
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
        std::cerr << "nsubmitted: " << submitted << std::endl;
    }
#if 1
    while (!std::all_of(stats.begin(), stats.end(), [](Statistic &s) { return s.done; })) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
    }
#endif
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