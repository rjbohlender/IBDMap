//
// Created by Bohlender,Ryan James on 9/4/19.
//

#ifndef CARVAIBD_PARSER_HPP
#define CARVAIBD_PARSER_HPP

#include <string>
#include <armadillo>
#include <utility>

#include "split.hpp"
#include "isgzipped.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/**
 * @brief Class to handle converting individual indices into row index
 *
 * Imposes an order on the sample indices, and maps individual id pairs to row indices.
 */
struct Indexer {
  // Class counts and categories
  arma::uword case_count;
  arma::uword cont_count;
  std::vector<int> phenotypes;
  std::map<int, int> ordered_positions;

  // Pairs
  arma::uword case_case;
  arma::uword case_cont;

  Indexer() = default;
  Indexer(arma::uword case_count_, arma::uword cont_count_, std::vector<int> phenotypes_)
  : case_count(case_count_), cont_count(cont_count_), phenotypes(phenotypes_) {
    case_case = case_count * (case_count - 1) / 2.;
    case_cont = case_count * cont_count;
    // TODO fully construct the ordered values
    int case_idx = 0;
    int cont_idx = 0;
    for(int i = 0; i < phenotypes.size(); ++i) {
      switch(phenotypes[i]) {
      case 1:
        ordered_positions[i] = case_idx;
        case_idx++;
        break;
      case 0:
		ordered_positions[i] = cont_idx;
		cont_idx++;
		break;
      default:
        throw(std::runtime_error("Invalid phenotype value."));
      }
    }
    assert(case_idx == case_count && cont_idx == cont_count);
#ifndef NDEBUG
	// Test indexer translate code
	std::vector<std::string> test_pack;
	std::vector<int> case_indices = {0, 1, 4, 6, 8, 10};
	std::vector<int> cont_indices = {2, 3, 5, 7, 9, 11};
	for (auto i : case_indices) {
	  for (auto j : case_indices) {
		if (j <= i) {
		  continue;
		}
		std::stringstream ss;
		ss << ordered_positions[i] << "," << ordered_positions[j];
		test_pack.push_back(ss.str());
	  }
	}
	for(auto i : case_indices) {
	  for(auto j : cont_indices) {
		std::stringstream ss;
		ss << ordered_positions[i] << "," << ordered_positions[j];
		test_pack.push_back(ss.str());
	  }
	}
	for (const auto &v : test_pack) {
	  std::cerr << v << " ";
	}
	std::cerr << std::endl;
	std::cerr << "Translate_test: 0,1: " << translate(0, 1) << " ; " << test_pack[translate(0, 1)] << std::endl;
	std::cerr << "Translate_test: 2,4: " << translate(2, 4) << " ; " <<  test_pack[translate(2, 4)] << std::endl;
	std::cerr << "Translate_test: 3,6: " << translate(3, 6) << " ; " <<  test_pack[translate(3, 6)] << std::endl;
	std::cerr << "Translate_test: 4,11: " << translate(4, 11) << " ; " <<  test_pack[translate(4, 11)] << std::endl;
#endif
  }

  arma::uword translate(arma::uword i, arma::uword j) {
    arma::uword row = 0;
    if (phenotypes[i] == 0) {
      if (phenotypes[j] == 0) {
		throw(std::runtime_error("ERROR: Trying to find a control:control pair."));
	  } else { // case-control or control-case either way
		row += case_case; // Shift to case-control pairs
		row += ordered_positions[j] * cont_count + ordered_positions[i]; // Every case pairs with every control, so we move for all of the former pairs, to the current case, and the current control;
	  }
	}
    if (phenotypes[j] == 0) { // Case-Control pair
      row += case_case; // Shift to case-control pairs
      row += ordered_positions[i] * cont_count + ordered_positions[j]; // Every case pairs with every control, so we move for all of the former pairs, to the current case, and the current control;
    } else { // Case-Case pair
      // No offset for case-case pairs, start from 0
      // Every subsequent case has fewer subsequent pairs
      int case1 = ordered_positions[i];
      int case2 = ordered_positions[j];
      for(int i = 0; i < case1; i++) {
        row += case_count - i; // Move to the current case
      }
      row += std::abs(case1 - case2);
    }
    return row;
  }
};

/**
 * @brief Basic POD type to hold breakpoint information
 */
struct Breakpoint {
  std::string breakpoint;
  std::vector<double> segment_lengths;
  std::vector<std::pair<unsigned long, unsigned long>> ibd_pairs;
};

/**
 * @brief Unified IBD data format parser
 *
 * Discussed format:
 * # Header lines prefixed with hash
 * <A> <0.5;1,2> <0.6:3,4>
 * <B> <cM distance:pair sharing> <cM distance: pair sharing>
 *
 * So each line is a breakpoint identifier, eg. chromosome:position, followed by white space separated groups
 * of data that are shared at that breakpoint. Each breakpoint, the list of shared
 * segments should change, by definition of a breakpoint. A point of transition
 * where some haplotypes end and others begin.
 *
 * Data may be gzipped or uncompressed. Gzip detection will be done via the first two bytes.
 *
 * TODO -- Clarify internal representation for analysis.
 * 		-- What sort of analysis do we want to do?
 * 		-- What sort of extensions to the initial statistics do we want to be able to do?
 *
 * 	Should need a case-control status vector, which will be permuted.
 * 	Will need to be able to quickly calculate the test statistic, which should be:
 * 	\[
 * 		(cc - \bar{cc}) - (c/c - \bar{c/c})
 * 	\]
 * 	Where cc is the rate of sharing at the locus for case-case pairs, and c/c is case-control pairs.
 *
 * 	We should be able to drop the genome-wide average $\bar{cc}$ and $\bar{c/c}$ during
 * 	permutation as, in the long run, they are equal to each other in expectation under
 * 	the null.
 *
 * 	A 2D matrix, with rows = sample pairs, columns = breakpoints, should be enough.
 * 	This matrix should be very sparse. Most pairs will not be IBD at every locus.
 * 	Size $\approx (n(n-1)/2 + n * m) \times c$ where n, m are the number of
 * 	cases and controls, and c is the number of SNPs. We have to use sparse matrices
 * 	for this, because a dense matrix is far too large. Assuming n = 100, m = 100, c = 1e6,
 * 	we have a 396GB dense matrix.
 *
 * 	Example from real data:
 * 	n = 8312
 * 	m = 29861
 *	c = 10261
 *	(8312 * (8311) / 2 + 8312 * 29861) * 10261 * 8 = 23209983709024.0 bytes
 *	== 23209.98 GB Dense matrix
 *	At ~0.01% fill rate - 2GB usage?
 *
 * 	We also need to order the sample pairs, so that we can quickly index to case-case and case-control
 * 	groups. The Indexer class will take care of this.
 *
 * 	If we do need to calculate the genome-wide average, then we could parallelize the work
 * 	across the matrix, with each thread handling a chunk of the data, then returning all the
 * 	results and merging them into the average.
 *
 * 	An alternative data format?
 * 	For each breakpoint -- The number of IBD pairs.
 * 	The total number of pairs is constant.
 * 	The statistic is based on the rate of sharing at the locus relative to the rate of sharing globally.
 * 	But how do we permute? How do we determine the case-control status of the individuals?
 *
 * 	TODO We have a problem with permutation of the test statistic
 * 	The statistic doesn't depend on control-control pairs. Under permutation, many of our original pairs will be
 * 	control control pairs. Is there a problem with pairs moving in and out of the sample space?
 */
template<class StringT>
class Parser {
  std::vector<std::string> samples;
  std::vector<int> phenotypes;
  Indexer indexer;
  std::map<int, std::string> sample_idx_map;
  arma::sp_mat data;
  std::vector<Breakpoint> breakpoints;

  void parse_input(std::istream &is) {
	std::string line;
	unsigned long varno = 0;
	while(std::getline(is, line)) {
	  RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	  if (boost::starts_with(line, "#")) {
	    if (splitter[1] == "breakpoints:") { // Allows me to set the size of the matrix without reading the whole file
	      data.set_size(indexer.case_case + indexer.case_cont, std::stoul(splitter[2]));
	      continue;
	    } else {
		  // Header
		  int index = std::stoi(splitter[1]);
		  std::string sample = splitter[2];
		  sample_idx_map[index] = sample;
		  continue;
		}
	  }
	  if (splitter.size() > 1) { // Possibility of empty entries (ending breakpoints)
	    for (int i = 1; i < splitter.size(); i++) {
	      RJBUtil::Splitter<std::string_view> entry(splitter[i], ";");
	      RJBUtil::Splitter<std::string_view> pairs(entry[1], ",");

	      unsigned long pair1 = std::stoul(pairs[0]);
		  unsigned long pair2 = std::stoul(pairs[1]);

		  if (i == 1) {
			Breakpoint bp {
				splitter[0],
				{std::stod(entry[0])},
				{std::make_pair(pair1, pair2)}
			};
			breakpoints.push_back(bp);
		  } else {
	        breakpoints.back().segment_lengths.push_back(std::stod(entry[0]));
			breakpoints.back().ibd_pairs.push_back(std::make_pair(pair1, pair2));
	      }

		  data(indexer.translate(pair1, pair2), varno) = 1;
	    }
	  } else {
	    Breakpoint bp {
	      splitter[0],
		  {},
		  {}
	    };
	    breakpoints.push_back(bp);
	  }
	  varno++;
	}
  }

  void parse_ped(std::istream &is) {
    // Format columns
    const int fid = 0;
    const int iid = 1;
    const int pid = 2;
    const int mid = 3;
    const int sex = 4;
    const int aff = 5;

	std::string line;
	while(std::getline(is, line)) {
	  if(boost::starts_with(line, "#")) { // Skip possible header
	    continue;
	  }
	  RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	  if (splitter.size() < 6) { // Ensure that we have enough elements to parse
	    throw(std::runtime_error("ERROR: Incorrectly formatted .ped file."));
	  }
	  samples.push_back(splitter[iid]);
	  phenotypes.push_back(std::stoi(splitter[aff]) - 1); // 2 is case, 1 is control, 0 is unknown
	}
	int case_count = std::accumulate(phenotypes.begin(), phenotypes.end(), 0);
	int control_count = phenotypes.size() - case_count;

	indexer = Indexer(case_count, control_count, phenotypes);
  }
public:
  /**
   * @brief Constructor for the Parser -- Handles all input parsing.
   * @param ipath Unified IBD format input file path
   * @param ppath Ped path
   */
  Parser(StringT ipath, StringT ppath) {
	IsGzipped<StringT> is_gzipped;
	std::ifstream ifs;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> streambuf;
	if (is_gzipped(ipath)) {
	  ifs.open(ipath, std::ios_base::in | std::ios_base::binary);
	  streambuf.push(boost::iostreams::gzip_decompressor());
	  streambuf.push(ifs);
	} else {
	  ifs.open(ipath, std::ios_base::in);
	  streambuf.push(ifs);
	}
	std::istream is(&streambuf);
	std::ifstream pifs(ppath);

	parse_ped(pifs);
	parse_input(is);
#ifndef NDEBUG
	std::cerr << data;
	std::cerr << data.rows(0, indexer.case_case - 1);
	std::cerr << data.rows(indexer.case_case, data.n_rows - 1);
	std::cerr << arma::mean(arma::mean(data.rows(0, indexer.case_case - 1))) << std::endl;
	std::cerr << arma::mean(arma::mean(data.rows(indexer.case_case, data.n_rows - 1))) << std::endl;
#endif
  }
};

#endif //CARVAIBD_PARSER_HPP
