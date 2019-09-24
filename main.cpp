#include <iostream>
#include <boost/program_options.hpp>
#include "src/parser.hpp"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  // C++ IO only
  std::ios_base::sync_with_stdio(false);

  po::options_description desc;

  desc.add_options()
		  ("input,i",
		   po::value<std::string>()->required(),
		   "Unified IBD input file.")
		  ("ped,p",
		   po::value<std::string>()->required(),
		   ".ped format file describing case-control status of all samples.")
		  ("help,h", "Display this message.");

  po::variables_map vm;
  try {
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	  std::cerr << desc << std::endl;
	  return 1;
	}
  } catch(po::error &pe) {
    std::cerr << pe.what() << std::endl;
    std::cerr << desc << std::endl;
    return 1;
  }
  Parser parser(vm["input"].as<std::string>(), vm["ped"].as<std::string>());
}