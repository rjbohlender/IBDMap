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
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/optional.hpp>

template<class StringT>
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
  void count_breakpoints(std::istream &is) {
	nbreakpoints = 0;
	std::string line;
	while (std::getline(is, line)) {
	  if (boost::starts_with(line, "#"))
		continue;
	  nbreakpoints++;
	}
	if (params.verbose) {
	  std::cerr << "total breakpoints: " << nbreakpoints << std::endl;
	}
  }

  void parse_input(std::istream &is) {
	std::string line;
	unsigned long submitted = 0;
	double cur_dist = 0;
	double last_dist = 0;

	// Initialize ThreadPool
	ThreadPool<Statistic> threadpool(params);

	// Maintain a single vector that we just update with each line
	arma::sp_vec last(samples.size() * (samples.size() - 1) / 2.);
	arma::sp_vec data(samples.size() * (samples.size() - 1) / 2.);

	// 3 pointers for each statistic class
	auto index_ptr = std::shared_ptr<std::vector<Indexer>>(&indexer);
	auto sample_ptr = std::shared_ptr<std::vector<std::string>>(&samples);
	auto pheno_ptr = std::shared_ptr<std::vector<std::vector<int>>>(&phenotypes);

	// Generate permutations if we have covariates
	Permute permute(params.seed);
	boost::optional<std::shared_ptr<std::vector<std::vector<int32_t>>>> bo_perms;
	if (covariates && phenotypes.size() == 1) {
	  auto permutation_ptr = std::make_shared<std::vector<std::vector<int32_t>>>();

	  arma::vec Y = arma::conv_to<arma::vec>::from(phenotypes[0]);
	  Binomial link;
	  GLM<Binomial> fit(*covariates, Y, link);
	  arma::vec odds = fit.mu_ / (1 - fit.mu_);

	  std::string lr_path = params.output_path.empty() ? "lr.txt" : params.output_path + ".lr.txt";

	  std::ofstream lr_out(lr_path);
	  lr_out << "Sample\tProb\tOdds" << std::endl;
	  for (int i = 0; i < odds.n_elem; i++) {
		lr_out << cov_samples[i] << "\t" << fit.mu_(i) << "\t" << odds(i) << std::endl;
	  }

	  if(params.verbose) {
		std::cerr << "LR Output: max pr: " << arma::max(fit.mu_) << " min pr: " << arma::min(fit.mu_) << std::endl;
		std::cerr << "LR Output: max odds: " << arma::max(odds) << " min odds: " << arma::min(odds) << std::endl;
	  }
	  lr_out.close();

	  permute.get_permutations(permutation_ptr, odds, indexer[0].case_count, params.nperms, params.nthreads - 2);
	  bo_perms = permutation_ptr;
	}

	arma::wall_clock timer;
	long lineno = -1;
	while (std::getline(is, line)) {
	  lineno++;
	  if (lineno == 0) { // Skip the header
		continue;
	  }
	  int cscs_cnt = 0;
	  int cscn_cnt = 0;

	  RJBUtil::Splitter<std::string_view> splitter(line, "\t", true);

	  std::string chrom;

	  if (!boost::starts_with(splitter[0], "chr")) {
		std::stringstream chrss;
		chrss << "chr" << splitter[0];
		chrom = chrss.str();
	  }

	  int pos = std::stoi(splitter[1]);

	  RJBUtil::Splitter<std::string_view> additions(splitter[splitter.size() - 2], " ");
	  RJBUtil::Splitter<std::string_view> deletions(splitter[splitter.size() - 1], " ");
	  // Initialize breakpoint
	  Breakpoint bp{};
	  bp.breakpoint = std::make_pair(chrom, splitter[1]);

	  for (auto &entry : additions) {
		if (entry == "NA") {
		  break;
		}
		RJBUtil::Splitter<std::string_view> vals(entry, ":");
		RJBUtil::Splitter<std::string> pairs(vals[1], "-");

		arma::sword row_idx = indexer[0].translate(pairs[0], pairs[1]);
		if (row_idx < 0) {
		  continue;
		} else {
		  try {
			bp.segment_lengths.push_back(std::stod(vals[0]));
		  } catch(std::invalid_argument &e) {
		    bp.segment_lengths.push_back(std::nan("1"));
		  }
		  bp.ibd_pairs.emplace_back(std::make_pair(pairs[0], pairs[1]));
		}
		data(row_idx) += 1;
	  }
	  for (auto &entry : deletions) {
		if (entry == "NA") {
		  break;
		}
		RJBUtil::Splitter<std::string_view> vals(entry, ":");
		RJBUtil::Splitter<std::string> pairs(vals[1], "-");

		arma::sword row_idx = indexer[0].translate(pairs[0], pairs[1]);
		if (row_idx < 0) {
		  continue;
		} else {
		  try {
			bp.segment_lengths.push_back(std::stod(vals[0]));
		  } catch(std::invalid_argument &e) {
			bp.segment_lengths.push_back(std::nan("1"));
		  }
		  bp.ibd_pairs.emplace_back(std::make_pair(pairs[0], pairs[1]));
		}
		data(row_idx) -= 1;
	  }

	  // Must follow data update -- Data format is just change from previous breakpoint so skipping update is impossible
	  std::pair<std::pair<int, double>, std::pair<int, double>> nearest = gmap.find_nearest(chrom, pos);
	  if (nearest.first.first != nearest.second.first) {
		cur_dist = (pos - nearest.first.first) * (nearest.second.second - nearest.first.second)
			/ (nearest.second.first - nearest.first.first) + nearest.first.second;
	  } else {
		cur_dist = nearest.first.second;
	  }

	  if (cur_dist - last_dist < params.min_dist) {
		continue;
	  }

	  auto cor = [](arma::sp_vec X, arma::sp_vec Y) -> double {
		double N = static_cast<double>(X.n_elem);
		double mxy = arma::as_scalar(arma::mean(X % Y)); // as_scalar requrired for compatibility with older versions
		double mx = arma::mean(X);
		double my = arma::mean(Y);
		double sdx = std::sqrt(arma::var(X));
		double sdy = std::sqrt(arma::var(Y));

		double r = N / (N - 1) * (mxy - mx * my) / (sdx * sdy);
		return std::pow(r, 2);
	  };

	  if (params.rsquared) {
		double r2 = cor(data, last);
		if (params.verbose) {
		  std::cerr << "r2: " << r2 << std::endl;
		}
		if (r2 > *params.rsquared) {
		  continue;
		}
	  }
	  if (params.verbose) {
		std::cerr << "dist: " << cur_dist - last_dist << std::endl;
		std::cerr << "fill rate: " << arma::accu(data) / data.n_elem << std::endl;
	  }
	  last_dist = cur_dist;
	  last = data;

	  Statistic stat(data,
					 bp,
					 index_ptr,
					 sample_ptr,
					 pheno_ptr,
					 reporter,
					 params,
					 groups,
					 bo_perms);
	  threadpool.submit(std::move(stat));
	  submitted++;
	  // for(auto it = stats.begin(); it != stats.end(); it++) { // Cleanup as we go
	// 	if ((*it).done) {
	// 	  it = stats.erase(it);
	// 	  std::cerr << "deleted\n";
	// 	}
	  // }
	}
	while (!threadpool.done) {
	  // for(auto it = stats.begin(); it != stats.end(); it++) { // Cleanup as we go
	// 	if ((*it).done) {
	// 	  it = stats.erase(it);
	// 	}
	  // }
	  std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
	}
  }

  void parse_pheno(std::istream &is) {
	int iid = 0;
	int phe = 1;
	std::string line;
	arma::uword lineno = 0;
	std::map<std::vector<bool>, std::vector<arma::uword>> fill_patterns;

	bool notified = false;

	while (std::getline(is, line)) {
	  if (boost::starts_with(line, "#") || lineno == 0) { // Skip the header
		lineno++;
		continue;
	  }
	  RJBUtil::Splitter<std::string_view> splitter(line, " \t");
	  if (splitter.size() > 2 && !notified) {
		std::cerr << "Multiple phenotypes provided. Covariates will be ignored." << std::endl;
		notified = true;
	  }

	  if (!splitter.empty() && skip.find(splitter[0]) != skip.end()) {
		// Skip samples with missing cov values; don't increment lineno because we're treating them as if they don't exist
		continue;
	  }

	  samples.push_back(splitter[iid]);
	  std::vector<bool> pattern;
	  for (int i = 1; i < splitter.size(); i++) {
		if (phenotypes.size() < i) {
		  phenotypes.emplace_back();
		}
		if (splitter[i] == "NA") {
		  pattern.push_back(false);
		  phenotypes[i - 1].push_back(-1);
		} else {
		  pattern.push_back(true);
		  phenotypes[i - 1].push_back(std::stoi(splitter[i]));
		  if (phenotypes[i - 1].back() != 0 && phenotypes[i - 1].back() != 1) {
		    if(params.verbose) {
			  std::cerr << splitter[0] << " " << splitter[1] << std::endl;
		    }
		  }
		  if (params.swap) { // Swap case-control status
			switch (phenotypes[i - 1].back()) {
			case 1:phenotypes[i - 1].back() = 0;
			  break;
			case 0:phenotypes[i - 1].back() = 1;
			  break;
			default:break;
			}
		  }
		}
	  }
	  if (fill_patterns.count(pattern) == 0) {
		fill_patterns.emplace(std::make_pair(pattern, std::vector<arma::uword>({lineno - 1})));
	  } else {
		fill_patterns[pattern].push_back(lineno - 1);
	  }
	  lineno++;
	}
	for (const auto &v : phenotypes) {
	  int case_count = 0;
	  int control_count = 0;
	  for (const auto &p : v) {
		switch (p) {
		case 1:case_count++;
		  break;
		case 0:control_count++;
		  break;
		default:break;
		}
	  }

	  indexer.emplace_back(Indexer(case_count, control_count, samples, phenotypes.back()));
	}
	if (fill_patterns.size() > 1) {
	  groups = std::vector<std::vector<arma::uword>>();
	  for (const auto &v : fill_patterns) {
		groups->push_back(v.second);
	  }
	  if(params.verbose) {
		std::cerr << "Groups: " << groups->size() << std::endl;
		std::cerr << "Group sizes: ";
	  }
	  for (const auto &v : fill_patterns) {
	    if (params.verbose) {
		  std::cerr << v.second.size() << " ";
	    }
	  }
	  if(params.verbose) {
		std::cerr << std::endl;
	  }
	}
	if (covariates) {
	  int i = 0;
	  arma::uvec idx(samples.size());
	  for (const auto &v : samples) {
		if (cov_samples[i] != v) {
		  auto swap_val = std::find(cov_samples.begin(), cov_samples.end(), v);
		  if (swap_val == cov_samples.end())
			throw (std::runtime_error("ERROR: Sample not present in covariate file."));
		  int j = std::distance(cov_samples.begin(), swap_val);
		  idx(i) = j;
		  std::swap(cov_samples[i], cov_samples[j]);
		} else {
		  idx(i) = i;
		}
		i++;
	  }
	  covariates = (*covariates).rows(idx);
	}
  }

  void parse_cov(std::istream &is) {
	std::string line;
	unsigned long lineno = 0;
	int nfields = 0;

	std::map<std::string, std::vector<double>> data;
	std::vector<std::vector<std::string>> unconvertible;

	while (std::getline(is, line)) {
	  RJBUtil::Splitter<std::string> splitter(line, " \t");
	  if (lineno == 0) { // Skip the header
		if (splitter.size() > 0) {
		  nfields = splitter.size() - 1; // Not counting the sample field so we don't have to subtract all the time.
		  unconvertible.resize(nfields);
		  lineno++;
		} else {
		  throw (std::runtime_error("Header line of covariate file is empty. Please include a header line."));
		}
		continue;
	  }

	  if (splitter.empty()) {
		continue;
	  }
	  std::string sampleid = splitter[0];
	  if (splitter.size() < nfields + 1) { // Skip samples where we have a missing column or two
		skip.emplace(sampleid);
		continue;
	  }

	  cov_samples.push_back(sampleid);
	  data[sampleid] = std::vector<double>(nfields, 0);

	  for (int i = 0; i < nfields; i++) {
		try {
		  data[sampleid][i] = std::stod(splitter[i + 1]);
		} catch (...) {
		  unconvertible[i].push_back(splitter[i + 1]);
		}
	  }
	  lineno++;
	}

	// Handle unconvertible fields by treating them as factors with levels -- convert to dummy variables
	int fieldno = 0;
	int offset = 0;
	for (const auto &field : unconvertible) {
	  if (!field.empty()) {
		std::set<std::string> unique(field.begin(), field.end());
		std::map<std::string, int> levels;

		if (params.verbose) {
		  std::cerr << "In reading covariates, could not convert column " << fieldno + 1 << " to double." << std::endl;
		  std::cerr << "Levels: ";
		}

		for (auto it = unique.begin(); it != unique.end(); it++) {
		  levels.emplace(std::make_pair(*it, std::distance(unique.begin(), it)));
		  if (params.verbose) {
			std::cerr << *it << " : " << std::distance(unique.begin(), it) << " ";
		  }
		}
		if (params.verbose) {
		  std::cerr << std::endl;
		}

		int sampleno = 0;
		int nlevels = levels.size() - 1;
		for (const auto &v : field) { // Convert to dummy variable
		  for (int j = 0; j < nlevels; j++) {
			if (j == 0) {
			  if (j == levels[v]) {
				data[cov_samples[sampleno]][fieldno + offset] = 1.0;
			  } else {
				data[cov_samples[sampleno]][fieldno + offset] = 0.0;
			  }
			} else {
			  if (j == levels[v]) {
				data[cov_samples[sampleno]].insert(data[cov_samples[sampleno]].begin() + fieldno + offset + j, 1.0);
			  } else {
				data[cov_samples[sampleno]].insert(data[cov_samples[sampleno]].begin() + fieldno + offset + j, 0.0);
			  }
			}
		  }
		  sampleno++;
		}
		offset += nlevels - 1;
		nfields += nlevels - 1;
	  }
	  fieldno++;
	}

	arma::mat design(cov_samples.size(), nfields + 1);
	if(params.verbose) {
	  std::cerr << "Design.n_rows: " << design.n_rows << std::endl;
	  std::cerr << "Design.n_cols: " << design.n_cols << std::endl;
	}
	int i = 0;
	for (const auto &s : cov_samples) {
	  design(i, 0) = 1;
	  int j = 1;
	  for (const auto &v : data[s]) {
		design(i, j) = v;
		j++;
	  }
	  i++;
	}

	covariates = design;
  }

public:
  std::vector<std::string> samples; // Samples in input order
  std::vector<std::string> cov_samples; // Samples in covariate order
  std::set<std::string> skip; // Skip samples with missing cov values
  std::vector<std::vector<int>> phenotypes;
  std::vector<Indexer> indexer;
  unsigned long nbreakpoints;
  Parameters params;
  std::shared_ptr<Reporter> reporter;
  boost::optional<std::vector<std::vector<arma::uword>>> groups;
  boost::optional<arma::mat> covariates;
  GeneticMap gmap;

  /**
   * @brief Parser and data dispatcher
   * @param input_path Path to the unified IBD format file
   * @param pheno_path Path to the phenotype description file
   * @param cov_path Path to the covariates -- for covariate adjusted permutation
   * @param params_ Parameter injection for program options
   * @param reporter_ Global reporter
   * @param gmap_ Genetic map file for calculation of distances between breakpoints
   */
  Parser(StringT input_path,
		 StringT pheno_path,
		 boost::optional<StringT> cov_path,
		 Parameters params_,
		 std::shared_ptr<Reporter> reporter_,
		 GeneticMap &gmap_)
	  : nbreakpoints(0), params(std::move(params_)), reporter(std::move(reporter_)), gmap(std::move(gmap_)) {
	Source bp_source(input_path);
	std::istream bp_is(&(*bp_source.streambuf));
	std::ifstream pheno_ifs(pheno_path);

	if(params.verbose) {
	  std::cerr << "Counting breakpoints\n";
	}
	count_breakpoints(bp_is);

	Source input_source(input_path);
	std::istream input_s(&(*input_source.streambuf));

	if (cov_path) {
	  if(params.verbose) {
		std::cerr << "Parsing covariates\n";
	  }

	  Source cov_source(*cov_path);
	  std::istream cov_is(&(*cov_source.streambuf));

	  parse_cov(cov_is);
	}

	if(params.verbose) {
	  std::cerr << "Parsing phenotypes\n";
	}
	parse_pheno(pheno_ifs);

	for (const auto &idx : indexer) {
	  if(params.verbose) {
		std::cerr << "ncases: " << idx.case_count << " ncontrols: " << idx.cont_count << std::endl;
	  }
	}

	if(params.verbose) {
	  std::cerr << "Parsing data\n";
	}
	parse_input(input_s);
  }
};

#endif //CARVAIBD_PARSER_HPP
