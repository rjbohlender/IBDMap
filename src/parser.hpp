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
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/**
 * @brief Class to handle converting individual indices into row index
 *
 * Invertible mapping function between sample pairs and rows.
 */
struct Indexer {
    // Class counts and categories
    arma::uword case_count;
    arma::uword cont_count;
    std::vector<int> phenotypes;
    std::vector<std::string> samples;
    std::vector<arma::sword> transitions; // Transition points between pairing sets
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
            : case_count(case_count_), cont_count(cont_count_), phenotypes(std::move(phenotypes_)),
              samples(std::move(samples_)) {
        setup(case_count_, cont_count_);
#ifndef NDEBUG
     test_run();
#endif
    }

    void setup(arma::uword case_count_,
               arma::uword cont_count_) {
        IndexSort indexSort(samples);
        indexSort.sort_vector(samples);
        indexSort.sort_vector(phenotypes); // Both must be sorted
        case_case = case_count * (case_count - 1.) / 2.;
        case_cont = case_count * cont_count;
        cont_cont = cont_count * (cont_count - 1.) / 2.;

        int case_idx = 0;
        int cont_idx = 0;
        int current = 0;

        arma::sword start = 0;
        arma::sword k = samples.size() - 1;
        for (; k > 0; k--) {
            start += k;
            transitions.push_back(start);
        }

        for (const auto &s : samples) {
            switch (phenotypes[current]) {
                case 1:
                    ordered_positions[s] = current;
                    ordered_cc_positions[s] = case_idx;
                    case_idx++;
                    break;
                case 0:
                    ordered_positions[s] = current;
                    ordered_cc_positions[s] = cont_idx;
                    cont_idx++;
                    break;
                default:
                    throw (std::runtime_error("Invalid phenotype value."));
            }
            current++;
        }
    }

    arma::sword translate(std::string a, std::string b) {
        auto finda = std::lower_bound(samples.begin(), samples.end(), a);
        auto findb = std::lower_bound(samples.begin(), samples.end(), b);
        if (*finda != a || *findb != b) {
            return -1;
        }
        int i = std::distance(samples.begin(), finda);
        int j = std::distance(samples.begin(), findb);
        if (i == j) {
            throw (std::runtime_error("Searching for two of the same pair."));
        }
        if (i > j) {
            auto tmp = i;
            i = j;
            j = tmp;
        }

        arma::sword start = 0;
        arma::sword k = samples.size() - 1;
        for (; k > samples.size() - i - 1; k--)
            start += k;
        return start + j - i - 1;
    }

    std::pair<arma::sword, arma::sword> back_translate(arma::uword row) {
        arma::sword i = -1;
        arma::sword j = -1;

        auto bound = std::lower_bound(transitions.begin(), transitions.end(), row);

        i = std::distance(transitions.begin(), bound);
        if(row == *bound) {
            i += 1;
            j = i + 1;
        } else {
            j = i > 0 ? row - *(bound - 1) + i + 1 : row + 1;
        }

        return std::make_pair(i, j);
    }

    void test_run() {
        std::vector<std::string> test_samples {
            "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"
        };

        std::vector<int> test_phenotypes {
            1, 1, 1, 1, 1, 0, 0, 0, 0, 0
        };

        arma::sword N = 10;
        std::vector<int> test_transitions;
        int start = 0;
        for(int k = N - 1; k > 0; k--) {
            start += k;
            test_transitions.push_back(start);
        }

        auto test_translate = [&](std::string &a, std::string &b) {
            auto finda = std::lower_bound(test_samples.begin(), test_samples.end(), a);
            auto findb = std::lower_bound(test_samples.begin(), test_samples.end(), b);
            if (*finda != a || *findb != b) {
                return -1ll;
            }
            int i = std::distance(test_samples.begin(), finda);
            int j = std::distance(test_samples.begin(), findb);
            std::cerr << "i: " << i << " j: " << j << std::endl;
            if (i == j) {
                throw (std::runtime_error("Searching for two of the same pair."));
            }
            if (i > j) {
                auto tmp = i;
                i = j;
                j = tmp;
            }

            arma::sword start = 0;
            arma::sword k = test_samples.size() - 1;
            for (; k > test_samples.size() - i - 1; k--)
                start += k;
            return start + j - i - 1;
        };

        auto test_back_translate = [&](arma::sword row) {
            arma::sword i = -1;
            arma::sword j = -1;

            auto bound = std::lower_bound(test_transitions.begin(), test_transitions.end(), row);

            i = std::distance(test_transitions.begin(), bound);
            if(row == *bound) {
                i += 1;
                j = i + 1;
            } else {
                j = i > 0 ? row - *(bound - 1) + i + 1 : row + 1;
            }
            std::cerr << "bound: " << *bound << std::endl;

            return std::make_pair(i, j);
        };

        // Testing A:B; Correct answer: 0
        auto row = test_translate(test_samples[0], test_samples[1]);
        std::cerr << "Row returned: " << row << " Correct row: " << 0 << std::endl;
        row = test_translate(test_samples[0], test_samples[2]);
        std::cerr << "Row returned: " << row << " Correct row: " << 1 << std::endl;
        row = test_translate(test_samples[0], test_samples[3]);
        std::cerr << "Row returned: " << row << " Correct row: " << 2 << std::endl;
        row = test_translate(test_samples[1], test_samples[2]);
        std::cerr << "Row returned: " << row << " Correct row: " << 9 << std::endl;
        row = test_translate(test_samples[8], test_samples[9]);
        std::cerr << "Row returned: " << row << " Correct row: " << 44 << std::endl;

        auto p = test_back_translate(0);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "0,1" << std::endl;
        p = test_back_translate(44);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "8,9" << std::endl;
        p = test_back_translate(9);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "1,2" << std::endl;
        p = test_back_translate(10);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "1,3" << std::endl;
        p = test_back_translate(5);
        std::cerr << "Pair returned: " << p.first << "," << p.second << " Correct pair: " << "0,6" << std::endl;
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

        data.set_size(indexer[0].case_case + indexer[0].case_cont + indexer[0].cont_cont, nbreakpoints);

        while (std::getline(is, line)) {
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
            std::cerr << "splitter size: " << splitter.size() << std::endl;
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
                            std::cerr << "mismatch: " << indexer[0].samples[bt.second] << " " << pair1  << " " << pair2 << std::endl;
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
                    if (data(row_idx, varno) > 0) {
                        std::cerr << "already written\n";
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
            std::cerr << splitter[0] << " " << splitter[1] << " cscs_cnt: " << cscs_cnt << " cscn_cnt: " << cscn_cnt
                      << std::endl;
            varno++;
        }
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
            samples.push_back(splitter[iid]);
            for (int i = 1; i < splitter.size(); i++) {
                if (phenotypes.size() < i) {
                    phenotypes.emplace_back();
                }
                if (splitter[i] == "NA") {
                    phenotypes[i - 1].push_back(-1);
                } else {
                    phenotypes[i - 1].push_back(std::stoi(splitter[i]));
                }
            }
        }
        for (unsigned long k = 0; k < phenotypes.size(); k++) {
            int case_count = 0;
            int control_count = 0;
            for (const auto &v : phenotypes[k]) {
                switch (v) {
                    case 1:
                        case_count++;
                        break;
                    case 0:
                        control_count++;
                        break;
                    default:
                        break;
                }
            }

            indexer.emplace_back(Indexer(case_count, control_count, samples, phenotypes[k]));
        }
    }

public:
    std::vector<std::string> samples;
    std::vector<std::vector<int>> phenotypes;
    std::vector<Indexer> indexer;
    arma::sp_mat data;
    std::vector<Breakpoint> breakpoints;
    unsigned long nbreakpoints;

    /**
     * @brief Constructor for the Parser -- Handles all input parsing.
     * @param ipath Unified IBD format input file path
     * @param ppath Phenotype path
     */
    Parser(StringT ipath, StringT ppath) : nbreakpoints(0) {
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
