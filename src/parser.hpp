//
// Created by Bohlender,Ryan James on 9/4/19.
//

#ifndef CARVAIBD_PARSER_HPP
#define CARVAIBD_PARSER_HPP

#include <string>
#include <armadillo>
#include <utility>
#include <unordered_map>
#include "absl/container/flat_hash_map.h"

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
  std::vector<std::string> samples;
  absl::flat_hash_map<arma::uword, std::pair<arma::sword, arma::sword>> row_map;
  absl::flat_hash_map<std::string, int> ordered_positions;
  absl::flat_hash_map<std::string, int> ordered_cc_positions;
  // std::unordered_map<std::string, int> ordered_positions;
  // std::unordered_map<std::string, int> ordered_cc_positions;


  // Pairs
  arma::uword case_case;
  arma::uword case_cont;
  arma::uword cont_cont;

  Indexer() = default;
  Indexer(arma::uword case_count_,
		  arma::uword cont_count_,
		  std::vector<std::string> samples_,
		  std::vector<int> phenotypes_)
	  : case_count(case_count_), cont_count(cont_count_), phenotypes(std::move(phenotypes_)), samples(std::move(samples_)) {
	setup(case_count_, cont_count_);
  }

  void setup(arma::uword case_count_,
			 arma::uword cont_count_) {
	case_case = case_count * (case_count - 1.) / 2.;
	case_cont = case_count * cont_count;
	cont_cont = cont_count * (cont_count - 1.) / 2.;

	int case_idx = 0;
	int cont_idx = 0;
	int current = 0;
	for (const auto &s : samples) {
	  switch (phenotypes[current]) {
	  case 1:ordered_positions[s] = current;
	  	ordered_cc_positions[s] = case_idx;
		case_idx++;
		break;
	  case 0: ordered_positions[s] = current;
		ordered_cc_positions[s] = cont_idx;
		cont_idx++;
		break;
	  default:throw (std::runtime_error("Invalid phenotype value."));
	  }
	  current++;
	}
#ifndef LOWMEM
	row_map.reserve(case_case + case_cont + cont_cont);

	std::cerr << "Starting to build row_map.\n";
	arma::wall_clock timer;
	timer.tic();
	arma::uword row = 0;
	for(arma::uword i = 0; i < samples.size(); i++) {
	  for(arma::uword j = i + 1; j < samples.size(); j++) {
	    row_map[row] = std::make_pair<arma::sword, arma::sword>(ordered_positions[samples[i]], ordered_positions[samples[j]]);
	    row++;
	  }
	}
	std::cerr << timer.toc() << " seconds to build row_map.\n";
#endif
	assert(case_idx == case_count && cont_idx == cont_count);
#ifndef NDEBUG
	// Test indexer translate code
	std::vector<std::string> test_pack;
	std::map<std::string, int> test_pos;
	std::map<std::string, int> test_cc_pos;
	std::vector<int> test_pheno = {1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0};
	std::vector<std::string> case_indices = {"0", "1", "4", "6", "8", "10"};
	std::vector<std::string> cont_indices = {"2", "3", "5", "7", "9", "11"};
	int test_case_case = case_indices.size() * (case_indices.size() - 1) / 2.;
	int test_case_cont = case_indices.size() * cont_indices.size();
	int test_cont_cont = cont_indices.size() * (cont_indices.size() - 1) / 2.;
	int test_case_count = 6;
	int test_cont_count = 6;

	// Fill test_pos
	case_idx = 0;
	cont_idx = 0;
	current = 0;
	for (const auto &phen : test_pheno) {
	  switch(phen) {
	  case 1: test_pos[case_indices[case_idx]] = current;
	  test_cc_pos[case_indices[case_idx]] = case_idx;
	  case_idx++;
	  break;
	  case 0: test_pos[cont_indices[cont_idx]] = current;
	  test_cc_pos[cont_indices[cont_idx]] = cont_idx;
	  cont_idx++;
	  break;
	  default: break;
	  }
	  current++;
	}

	// Test function
	auto test_translate = [&](const char* a, const char* b) {
	  int i = test_pos[a];
	  int j = test_pos[b];
	  arma::uword row = 0;
	  if (test_pheno[i] == 0) {
		if (test_pheno[j] == 0) {
		  row += test_case_case + test_case_cont;
		  int cont1 = test_cc_pos[a];
		  int cont2 = test_cc_pos[b] - 1;
		  for (int k = 1; k <= cont1; k++) {
			row += test_cont_count - k; // Move to the current case
		  }
		  row += std::abs(cont1 - cont2);
		} else { // control-case
		  row += test_case_case; // Shift to case-control pairs
		  row += test_cc_pos[b] * test_cont_count
			  + test_cc_pos[a]; // Every case pairs with every control, so we move for all of the former pairs, to the current case, and the current control;
		}
	  } else {
		if (test_pheno[j] == 0) { // Case-Control pair
		  row += test_case_case; // Shift to case-control pairs
		  row += test_cc_pos[a] * test_cont_count
			  + test_cc_pos[b]; // Every case pairs with every control, so we move for all of the former pairs, to the current case, and the current control;
		} else { // Case-Case pair
		  // No offset for case-case pairs, start from 0
		  // Every subsequent case has fewer subsequent pairs
		  int case1 = test_cc_pos[a];
		  int case2 = test_cc_pos[b] - 1;
		  for (int k = 1; k <= case1; k++) {
			row += test_case_count - k; // Move to the current case
		  }
		  row += std::abs(case1 - case2);
		}
	  }
	  return row;
	};

	for (int i = 0; i < case_indices.size(); i++) {
	  for (int j = 0; j < case_indices.size(); j++) {
		if (j <= i) {
		  continue;
		}
		std::stringstream ss;
		ss << test_pos[case_indices[i]] << "," << test_pos[case_indices[j]];
		test_pack.push_back(ss.str());
	  }
	}
	for (int i = 0; i < case_indices.size(); i++) {
	  for (int j = 0; j < cont_indices.size(); j++) {
		std::stringstream ss;
		ss << test_pos[case_indices[i]] << "," << test_pos[cont_indices[j]];
		test_pack.push_back(ss.str());
	  }
	}
	for (int i = 0; i < cont_indices.size(); i++) {
	  for (int j = 0; j < cont_indices.size(); j++) {
		if (j <= i) {
		  continue;
		}
		std::stringstream ss;
		ss << test_pos[cont_indices[i]] << "," << test_pos[cont_indices[j]];
		test_pack.push_back(ss.str());
	  }
	}
	for (const auto &v : test_pack) {
	  std::cerr << v << " ";
	}
	std::cerr << std::endl;
	std::cerr << "Test pack size: " << test_pack.size() << std::endl;
	std::cerr << "Translate_test: 0,1: " << test_translate("0", "1") << " ; " << test_pack[test_translate("0", "1")] << std::endl;
	std::cerr << "Translate_test: 2,4: " << test_translate("2", "4") << " ; " << test_pack[test_translate("2", "4")] << std::endl;
	std::cerr << "Translate_test: 3,6: " << test_translate("3", "6") << " ; " << test_pack[test_translate("3", "6")] << std::endl;
	std::cerr << "Translate_test: 4,11: " << test_translate("4", "11") << " ; " << test_pack[test_translate("4", "11")] << std::endl;
	std::cerr << "Translate_test: 9,11: " << test_translate("9", "11") << " ; " << test_pack[test_translate("9", "11")]
			  << std::endl;
#endif
  }

  arma::sword translate(std::string a, std::string b) {
    if (ordered_positions.count(a) < 1 || ordered_positions.count(b) < 1) {
      return -1;
    }

	int i = ordered_positions[a];
	int j = ordered_positions[b];
	arma::uword row = 0;
	if (phenotypes[i] == 0) {
	  if (phenotypes[j] == 0) {
		row += case_case + case_cont;
		int cont1 = ordered_cc_positions[a];
		int cont2 = ordered_cc_positions[b] - 1;
		if (cont1 < cont2) {
		  for (int k = 1; k <= cont1; k++) {
			row += cont_count - k; // Move to the current case
		  }
		} else {
		  for (int k = 1; k <= cont2; k++) {
			row += cont_count - k; // Move to the current case
		  }
		}
		row += std::abs(cont1 - cont2);
	  } else { // case-control or control-case either way
		row += case_case; // Shift to case-control pairs
		row += ordered_cc_positions[b] * cont_count
			+ ordered_cc_positions[a]; // Every case pairs with every control, so we move for all of the former pairs, to the current case, and the current control;
	  }
	} else {
	  if (phenotypes[j] == 0) { // Case-Control pair
		row += case_case; // Shift to case-control pairs
		row += ordered_cc_positions[a] * cont_count
			+ ordered_cc_positions[b]; // Every case pairs with every control, so we move for all of the former pairs, to the current case, and the current control;
	  } else { // Case-Case pair
		// No offset for case-case pairs, start from 0
		// Every subsequent case has fewer subsequent pairs
		int case1 = ordered_cc_positions[a];
		int case2 = ordered_cc_positions[b] - 1;
		if (case1 < case2) {
		  for (int k = 1; k <= case1; k++) {
			row += case_count - k; // Move to the current case
		  }
		} else {
		  for (int k = 1; k <= case2; k++) {
			row += case_count - k; // Move to the current case
		  }
		}
		row += std::abs(case1 - case2);
	  }
	}
	return row;
  }

  std::pair<arma::sword, arma::sword> back_translate(arma::uword row) {
    return row_map[row];
  }
  std::pair<arma::sword, arma::sword> back_translate_alt(arma::uword row) {
    if(row < case_case) {
      arma::sword case1 = 0, case2 = 0;
      arma::uword cur = 0;
      arma::sword k = 1;
      while(cur < row) {
        if (row - cur > case_count - k) {
          cur += case_count - k;
          k++;
          continue;
        }
        case1 = k;
        case2 = row - cur + k;
        if(case2 < 0) {
          std::cerr << "uhoh1\n";
        }
		if(ordered_positions[samples[case1]] != case1 || ordered_positions[samples[case2]] != case2) {
		  std::cerr << "ERROR: Failed to back translate.\n";
		}
		return std::make_pair(case1, case2);
      }
    } else if(row < case_case + case_cont) {
      arma::sword case1 = 0, cont1 = 0;
      arma::sword cur = case_case;
      arma::sword k = 0;

      while(cur < row) {
		if (row - cur > cont_count) {
		  cur += cont_count;
		  k++;
		  continue;
		}
		case1 = k;
		cont1 = row - cur;
		if(cont1 < 0) {
		  std::cerr << "uhoh2\n";
		}
		return std::make_pair(case1, cont1);
	  }
	} else if(row < case_case + case_cont + cont_cont) {
	  arma::sword cont1 = 0, cont2 = 0;
	  arma::uword cur = case_case + case_cont;
	  arma::sword k = 1;
	  while(cur < row) {
		if (row - cur > cont_count - k) {
		  cur += cont_count - k;
		  k++;
		  continue;
		}
		cont1 = k;
		cont2 = row - cur + k;
		if(cont2 < 0) {
		  std::cerr << "uhoh3\n";
		}
		return std::make_pair(cont1, cont2);
	  }
    } else {
      throw(std::runtime_error("ERROR: row value too large in back_translate_alt."));
    }
	throw(std::runtime_error("ERROR: back_translate_alt failed to find the pair."));
  }
};

/**
 * @brief Basic POD type to hold breakpoint information
 */
struct Breakpoint {
  std::pair<std::string, std::string> breakpoint;
  std::vector<double> segment_lengths;
  std::vector<std::pair<std::string, std::string>> ibd_pairs;
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
  std::map<int, std::string> sample_idx_map;

  struct Source {
	std::ifstream ifs;
	std::unique_ptr<boost::iostreams::filtering_streambuf<boost::iostreams::input>> streambuf;
	explicit Source(StringT path) {
	  IsGzipped<StringT> is_gzipped;
	  streambuf = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
	  if (is_gzipped(path)) {
		ifs.open(path, std::ios_base::in | std::ios_base::binary);
		(*streambuf).push(boost::iostreams::gzip_decompressor());
		(*streambuf).push(ifs);
	  } else {
		ifs.open(path, std::ios_base::in);
		(*streambuf).push(ifs);
	  }
	}
  };

  void count_breakpoints(std::istream &is) {
	nbreakpoints = 0;
	std::string line;
	while (std::getline(is, line)) {
	  if (boost::starts_with(line, "#"))
		continue;
	  nbreakpoints++;
	}
  }

  void parse_input(std::istream &is) {
	std::string line;
	unsigned long varno = 0;

	data.set_size(indexer.case_case + indexer.case_cont + indexer.cont_cont, nbreakpoints);

	while (std::getline(is, line)) {
	  RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	  if (boost::starts_with(line, "#")) {
		// Header
		int index = std::stoi(splitter[1]);
		std::string sample = splitter[2];
		sample_idx_map[index] = sample;
		continue;
	  }
	  if (splitter.size() > 2) { // Possibility of empty entries (ending breakpoints)
		for (unsigned long i = 2; i < splitter.size(); i++) {
		  RJBUtil::Splitter<std::string_view> entry(splitter[i], ":");
		  RJBUtil::Splitter<std::string_view> pairs(entry[1], "-");

		  std::string pair1(pairs[0]);
		  std::string pair2(pairs[1]);

		  arma::sword row_idx = indexer.translate(pair1, pair2);

		  if (i == 2) {
			Breakpoint bp{};
			if (row_idx < 0) {
			  bp.breakpoint = std::make_pair(splitter[0], splitter[1]);
			  breakpoints.push_back(bp);
			  continue;
			} else {
			  bp.breakpoint = std::make_pair(splitter[0], splitter[1]);
			  bp.segment_lengths.push_back(std::stod(entry[0]));
			  bp.ibd_pairs.emplace_back(std::make_pair(pair1, pair2));
			  breakpoints.push_back(bp);
			}
		  } else {
		    if (row_idx < 0) {
		      continue;
		    } else {
			  breakpoints.back().segment_lengths.push_back(std::stod(entry[0]));
			  breakpoints.back().ibd_pairs.push_back(std::make_pair(pair1, pair2));
		    }
		  }
		  data(row_idx, varno) = 1;
		}
	  } else {
		Breakpoint bp{
			std::make_pair(splitter[0], splitter[1]),
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
	while (std::getline(is, line)) {
	  if (boost::starts_with(line, "#")) { // Skip possible header
		continue;
	  }
	  RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	  if (splitter.size() < 6) { // Ensure that we have enough elements to parse
		throw (std::runtime_error("ERROR: Incorrectly formatted .ped file."));
	  }
	  samples.push_back(splitter[iid]);
	  phenotypes.push_back(std::stoi(splitter[aff]) - 1); // 2 is case, 1 is control, 0 is unknown
	}
	int case_count = std::accumulate(phenotypes.begin(), phenotypes.end(), 0);
	int control_count = phenotypes.size() - case_count;

	indexer = Indexer(case_count, control_count, phenotypes);
  }
  void parse_pheno(std::istream &is) {
	int iid = 0;
	int phe = 1;
	std::string line;
	int lineno = 0;
	while (std::getline(is, line)) {
	  if (boost::starts_with(line, "#") || lineno == 0) { // Skip the header
	    lineno++;
		continue;
	  }
	  RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	  if (splitter[1] == "NA")
	    continue;
	  samples.push_back(splitter[iid]);
	  phenotypes.push_back(std::stoi(splitter[phe]));
	}
	int case_count = std::accumulate(phenotypes.begin(), phenotypes.end(), 0);
	int control_count = phenotypes.size() - case_count;

	indexer = Indexer(case_count, control_count, samples, phenotypes);
  }

  void parse_avg(std::istream &is) {
    int iid1 = 0;
    int iid2 = 1;
    double avg = 0;
    std::string line;
    int lineno = 0;
    while(std::getline(is, line)) {
      if(boost::starts_with(line, "id1")) {
		continue;
	  }

      RJBUtil::Splitter<std::string> splitter(line, " \t");
      if(splitter.size() < 3) {
		throw std::runtime_error("ERROR: Incomplete line in pairwise avg file.\n");
	  }

      arma::sword row = indexer.translate(splitter[0], splitter[1]);

      if(row < 0) { // Pair not present in this phenotype set
        continue;
      }

      avg = std::stod(splitter[2]);
      row_pair_avg[row] = avg;
    }
  }
public:
  std::vector<std::string> samples;
  std::vector<int> phenotypes;
  Indexer indexer;
  arma::sp_mat data;
  std::vector<Breakpoint> breakpoints;
  unsigned long nbreakpoints;
  absl::flat_hash_map<arma::sword, double> row_pair_avg; // Maps the row of the pair to their genomewide avg
  /**
   * @brief Constructor for the Parser -- Handles all input parsing.
   * @param ipath Unified IBD format input file path
   * @param ppath Phenotype path
   */
  Parser(StringT ipath, StringT ppath, StringT apath) : nbreakpoints(0) {
	Source csource(ipath);
	std::istream cis(&(*csource.streambuf));
	std::ifstream pifs(ppath);

	count_breakpoints(cis);

	Source isource(ipath);
	std::istream is(&(*isource.streambuf));

	Source asource(apath);
	std::istream ais(&(*asource.streambuf));

	parse_pheno(pifs);
	parse_input(is);
	parse_avg(ais);

	std::cerr << "ncases: " << indexer.case_count << " ncontrols: " << indexer.cont_count << std::endl;
	std::cerr << "Pairs -- Case-Case: " << indexer.case_case << " Case-Control: " << indexer.case_cont << std::endl;
  }
};

#endif //CARVAIBD_PARSER_HPP
