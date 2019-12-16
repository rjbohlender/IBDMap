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
#include "IndexSort.hpp"
#include "parameters.hpp"
#include "reporter.hpp"
#include "threadpool.hpp"
#include "statistic.hpp"
#include "breakpoint.hpp"
#include "indexer.hpp"
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
  std::map<int, std::string> sample_idx_map;


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
    unsigned long submitted = 0;

    // Initialize ThreadPool
    ThreadPool<void> threadpool(params);
    std::vector<Statistic> stats;
    stats.reserve(nbreakpoints);

    while (std::getline(is, line)) {
      arma::sp_vec data(indexer[0].case_case + indexer[0].case_cont + indexer[0].cont_cont);

      int cscs_cnt = 0;
      int cscn_cnt = 0;

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

          arma::sword row_idx = indexer[0].translate(pair1, pair2);
          if (row_idx >= 0) {
            auto bt = indexer[0].back_translate(row_idx);
            if (indexer[0].samples[bt.first] != pair1 && indexer[0].samples[bt.first] != pair2) {
              std::cerr << "mismatch: " << indexer[0].samples[bt.first] << " " << pair1 << " " << pair2 << std::endl;
            }
            if (indexer[0].samples[bt.second] != pair1 && indexer[0].samples[bt.second] != pair2) {
              std::cerr << "mismatch: " << indexer[0].samples[bt.second] << " " << pair1 << " " << pair2 << std::endl;
            }
            if (phenotypes[0][bt.first] == 1) {
              if (phenotypes[0][bt.second] == 1) {
                cscs_cnt++;
              } else if (phenotypes[0][bt.second] == 0) {
                cscn_cnt++;
              }
            } else if (phenotypes[0][bt.first] == 0) {
              if (phenotypes[0][bt.second] == 1) {
                cscn_cnt++;
              }
            }
          }
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
          data(row_idx) = 1;
        }
      } else {
        Breakpoint bp{
            std::make_pair(splitter[0], splitter[1]),
            {},
            {}
        };
        breakpoints.push_back(bp);
      }

      Statistic stat(std::move(data),
                     breakpoints.back(),
                     indexer,
                     samples,
                     phenotypes,
                     reporter,
                     params);
        stats.emplace_back(std::move(stat));
        std::packaged_task<void()> f(std::bind(&Statistic::run, &stats.back()));
        while(threadpool.pending() >= 2 * params.nthreads) {
          std::this_thread::sleep_for(std::chrono::nanoseconds(1000000));
        }
        threadpool.submit(std::move(f));
        submitted++;
        std::cerr << "nsubmitted: " << submitted << std::endl;
      }
    while (!std::all_of(stats.begin(), stats.end(), [](Statistic &s) { return s.done; })) {
      std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
    }
  }

  void parse_pheno(std::istream &is) {
    int iid = 0;
    int phe = 1;
    std::string line;
    arma::uword lineno = 0;
    std::map<std::vector<bool>, std::vector<arma::uword>> fill_patterns;
    while (std::getline(is, line)) {
      if (boost::starts_with(line, "#") || lineno == 0) { // Skip the header
        lineno++;
        continue;
      }
      RJBUtil::Splitter<std::string_view> splitter(line, " \t");
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
        }
      }
      if (fill_patterns.count(pattern) == 0) {
        fill_patterns.emplace(std::make_pair(pattern, std::vector<arma::uword>({lineno})));
      } else {
        fill_patterns[pattern].push_back(lineno);
      }
      lineno++;
    }
    for (unsigned long k = 0; k < phenotypes.size(); k++) {
      int case_count = 0;
      int control_count = 0;
      for (const auto &v : phenotypes[k]) {
        switch (v) {
        case 1:case_count++;
          break;
        case 0:control_count++;
          break;
        default:break;
        }
      }

      indexer.emplace_back(Indexer(case_count, control_count, samples, phenotypes[k]));
    }
    if (fill_patterns.size() > 1) {
      for (const auto &v : fill_patterns) {
        groups->push_back(v.second);
      }
    }
  }

public:
  std::vector<std::string> samples;
  std::vector<std::vector<int>> phenotypes;
  std::vector<Indexer> indexer;
  std::vector<Breakpoint> breakpoints;
  unsigned long nbreakpoints;
  Parameters params;
  std::shared_ptr<Reporter> reporter;
  boost::optional<std::vector<std::vector<arma::uword>>> groups;

  /**
   * @brief Constructor for the Parser -- Handles all input parsing.
   * @param ipath Unified IBD format input file path
   * @param ppath Phenotype path
   */
  Parser(StringT ipath, StringT ppath, Parameters params_, std::shared_ptr<Reporter> reporter_)
      : nbreakpoints(0), params(params_), reporter(reporter_) {
    Source csource(ipath);
    std::istream cis(&(*csource.streambuf));
    std::ifstream pifs(ppath);

    std::cerr << "Counting breakpoints\n";
    count_breakpoints(cis);

    Source isource(ipath);
    std::istream is(&(*isource.streambuf));

    std::cerr << "Parsing phenotypes\n";
    parse_pheno(pifs);
    std::cerr << "Parsing data\n";
    parse_input(is);

    for (unsigned long k = 0; k < indexer.size(); k++) {
      std::cerr << "ncases: " << indexer[k].case_count << " ncontrols: " << indexer[k].cont_count << std::endl;
      std::cerr << "Pairs -- Case-Case: " << indexer[k].case_case << " Case-Control: " << indexer[k].case_cont
                << std::endl;
    }
  }
};

#endif //CARVAIBD_PARSER_HPP
