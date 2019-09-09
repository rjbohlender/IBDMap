//
// Created by Bohlender,Ryan James on 9/4/19.
//

#ifndef CARVAIBD_PARSER_HPP
#define CARVAIBD_PARSER_HPP

#include <string>
#include <armadillo>

/**
 * @brief Class to handle converting individual indices into row index
 *
 * Imposes an order on the sample indices, and maps individual id pairs to row indices.
 */
class Indexer {
  const arma::uword case_count;
  const arma::uword cont_count;
  arma::uword case_case;
  arma::uword case_cont;
public:
  Indexer(arma::uword case_count_, arma::uword cont_count_);

  arma::uword translate(arma::uword i, arma::uword j);
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
 * 	Size $\approx n^2m(n-1)/2 \times c$ where n, m are the number of
 * 	cases and controls, and c is the number of SNPs. We have to use sparse matrices
 * 	for this, because a dense matrix is far too large. Assuming n = 100, m = 100, c = 1e6,
 * 	we have a 396GB dense matrix.
 *
 * 	We also need to order the sample pairs, so that we can quickly index to case-case and case-control
 * 	groups.
 *
 * 	If we do need to calculate the genome-wide average, then we could parallelize the work
 * 	across the matrix, with each thread handling a chunk of the data, then returning all the
 * 	results and merging them into the average.
 */
class Parser {
  arma::sp_mat data;

  void parse(std::istream &is);
public:
  explicit Parser(const std::string &fpath);
  explicit Parser(std::string_view fpath);
};

#endif //CARVAIBD_PARSER_HPP
