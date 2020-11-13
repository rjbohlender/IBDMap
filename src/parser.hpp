//
// Created by Bohlender,Ryan James on 9/4/19.
//

#ifndef CARVAIBD_PARSER_HPP
#define CARVAIBD_PARSER_HPP

#define ARMA_DONT_USE_WRAPPER

#include <string>
#include <armadillo>
#include <utility>
#include <unordered_map>
#include <optional>
#include <set>

#include "split.hpp"
#include "isgzipped.hpp"
#include "indexsort.hpp"
#include "parameters.hpp"
#include "reporter.hpp"
#include "threadpool.hpp"
#include "statistic.hpp"
#include "breakpoint.hpp"
#include "indexer.hpp"
#include "permutation.hpp"
#include "../link/binomial.hpp"
#include "glm.hpp"
#include "geneticmap.hpp"
#include "info.hpp"
#include "math.hpp"
#include "source.hpp"
#include <boost/algorithm/string/predicate.hpp>

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
 */
class Parser {
  void count_breakpoints(std::istream &is);
  void parse_input(std::istream &is);
  void generate_cov_adj_perms(Permute &permute,
							  std::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> &o_perms);
  bool check_sample_list(const std::string &sample_pair);
  void update_data(arma::sp_vec &data,
				   std::map<std::string, int> &indices,
				   RJBUtil::Splitter<std::string_view> &changes,
				   Breakpoint &bp,
				   int value);
  void parse_pheno(std::istream &is);
  void parse_cov(std::istream &is);

public:
  std::shared_ptr<std::vector<std::string>> samples; // Samples in input order
  std::vector<std::string> cov_samples; // Samples in covariate order
  std::set<std::string> skip; // Skip samples with missing cov values
  std::vector<std::vector<int>> phenotypes;
  std::shared_ptr<std::vector<Indexer>> indexer;
  unsigned long nbreakpoints;
  Parameters params;
  std::shared_ptr<Reporter> reporter;
  std::optional<std::vector<std::vector<arma::uword>>> groups;
  std::optional<arma::mat> covariates;
  GeneticMap gmap;
  std::optional<Info> info;

  /**
   * @brief Parser and data dispatcher
   * @param input_path Path to the unified IBD format file
   * @param pheno_path Path to the phenotype description file
   * @param cov_path Path to the covariates -- for covariate adjusted permutation
   * @param params_ Parameter injection for program options
   * @param reporter_ Global reporter
   * @param gmap_ Genetic map file for calculation of distances between breakpoints
   */
  Parser(const std::string& input_path,
		 const std::string& pheno_path,
		 std::optional<std::string> cov_path,
		 Parameters params_,
		 std::shared_ptr<Reporter> reporter_,
		 GeneticMap &gmap_);
};

#endif //CARVAIBD_PARSER_HPP
