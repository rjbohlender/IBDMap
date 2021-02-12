//
// Created by Bohlender,Ryan James on 2/11/21.
//

#include <iostream>
#include <fmt/include/fmt/ostream.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "inputvalidator.hpp"
#include "split.hpp"

InputValidator::InputValidator(bool dash)
	: dash(dash) {}

void InputValidator::check_gmap(const std::string &line, size_t line_no) {
  std::string line_trim = boost::trim_right_copy(line);
  RJBUtil::Splitter<std::string> splitter(line_trim, " \t");
  if (splitter.size() != 3) {
	fmt::print(std::cerr, "Incorrect number of columns at line {} in GMAP input.", line_no);
	std::exit(-1);
  }
  if (line_no == 0) {
    if (splitter[0] != "pos") {
      fmt::print(std::cerr, "Incorrect GMAP header. Columns should be pos chr cM, separated by tabs. The header is expected to have this format.");
      std::exit(-1);
    }
    if (splitter[1] != "chr") {
	  fmt::print(std::cerr, "Incorrect GMAP header. Columns should be pos chr cM, separated by tabs. The header is expected to have this format.");
	  std::exit(-1);
	}
    if (splitter[2] != "cM") {
	  fmt::print(std::cerr, "Incorrect GMAP header. Columns should be pos chr cM, separated by tabs. The header is expected to have this format.");
	  std::exit(-1);
	}
  } else {
	try {
	  std::stoi(splitter[1]);
	} catch (std::exception &e) {
	  fmt::print(std::cerr, "Failed to parse position at line {} in GMAP input.", line_no);
	  std::exit(-1);
	}
	try {
	  std::stod(splitter[2]);
	} catch (std::exception &e) {
	  fmt::print(std::cerr, "Failed to parse CM at line {} in GMAP input.", line_no);
	  std::exit(-1);
	}
  }
}

void InputValidator::check_pheno(const std::string &line, size_t line_no) {
  RJBUtil::Splitter<std::string> splitter(line, " \t");
  if (line_no == 0) {
	if (splitter.size() < 2) {
	  fmt::print(std::cerr, "Phenotype file must have at least 2 columns. The sampleid and the phenotype.");
	  std::exit(-1);
	}
	pheno_column_count = splitter.size();
  }
  if (splitter.size() != pheno_column_count) {
	fmt::print(std::cerr, "Phenotype file has an incorrect number of columns at line {}.", line_no);
	std::exit(-1);
  }
}

void InputValidator::check_ibd(const std::string &line, size_t line_no) {
  RJBUtil::Splitter<std::string> splitter(line, "\t");
  if (dash) {
	if (boost::starts_with(line, "#")) {
	  // In header
	  return;
	}
	if (line_no == 0) {
	  // Header should look like: chr     pos     n.cluster       n.cluster.haplotype     n.cluster.pair  n.singleton.pair        cluster.add     cluster.del     singleton.add   singleton.del
	  // All tab separated
	  std::vector<std::string> columns{
		  "chr", "pos", "n.cluster", "n.cluster.haplotype", "n.cluster.pair", "n.singleton.pair",
		  "cluster.add", "cluster.del", "singleton.add", "singleton.del"
	  };

	  if (splitter.size() != columns.size()) {
		std::vector<std::string> missing;
		for (const auto &s : columns) {
		  if (std::find(splitter.begin(), splitter.end(), s) == splitter.end()) {
			missing.push_back(s);
		  }
		}

		fmt::print(std::cerr,
				   "IBD input header not properly formatted. Missing columns: {}",
				   fmt::join(missing.begin(), missing.end(), " "));
		std::exit(-1);
	  } else {
		ibd_column_count = splitter.size();
	  }
	} else {
	  if (splitter.size() != ibd_column_count) {
		fmt::print(std::cerr,
				   "Incorrect column count in IBD input at line {}. There are {} columns, but there should be {}.",
				   line_no,
				   splitter.size(),
				   ibd_column_count);
		std::exit(-1);
	  }
	}
  } else {
	if (line_no == 0) {
	  // Header should look like: chr     pos     segments        pairs   add     del
	  // All tab separated
	  std::vector<std::string> columns{
		  "chr", "pos", "segments", "pairs", "add", "del"
	  };
	  if (splitter.size() != columns.size()) {
		std::vector<std::string> missing;
		for (const auto &s : columns) {
		  if (std::find(splitter.begin(), splitter.end(), s) == splitter.end()) {
			missing.push_back(s);
		  }
		}

		fmt::print(std::cerr,
				   "IBD input header not properly formatted. Missing columns: {}",
				   fmt::join(missing.begin(), missing.end(), " "));
		std::exit(-1);
	  } else {
		ibd_column_count = splitter.size();
	  }
	} else {
	  if (splitter.size() != ibd_column_count) {
		fmt::print(std::cerr,
				   "Incorrect column count in IBD input at line {}. There are {} columns, but there should be {}.",
				   line_no,
				   splitter.size(),
				   ibd_column_count);
		std::exit(-1);
	  }
	}
  }
}

void InputValidator::check_info(const std::string &line, size_t line_no, const bool AF, const bool cM) {
  std::string line_trim = boost::trim_right_copy(line);
  RJBUtil::Splitter<std::string> splitter(line_trim, " \t");
  if (line_no == 0) {
    // Expected header: chr     clusterID       start   end     cM      count   freq
	std::vector<std::string> columns {
	  "chr", "clusterID", "start", "end", "count"
	};
	if (AF) {
	  columns.emplace_back("freq");
	}
	if (cM) {
	  columns.emplace_back("cM");
	}
	std::vector<std::string> missing;
	for (const auto &s : columns) {
	  if (std::find(splitter.begin(), splitter.end(), s) == splitter.end()) {
		missing.push_back(s);
	  }
	}

	if (missing.size() > 0) {
	  fmt::print(std::cerr,
				 "Info input header not properly formatted. Missing columns: {}",
				 fmt::join(missing.begin(), missing.end(), " "));
	  std::exit(-1);
	} else {
	  info_column_count = splitter.size();
	}
  } else {
    if (splitter.size() != info_column_count) {
      fmt::print(std::cerr, "Incorrect number of columns in info file at line {}.", line_no);
      std::exit(-1);
    }
  }
}
