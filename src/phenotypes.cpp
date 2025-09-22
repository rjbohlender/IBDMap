//
// Created by Bohlender,Ryan J on 10/8/21.
//

#include "phenotypes.hpp"
#include "../link/binomial.hpp"
#include "glm.hpp"
#include "permutation.hpp"
#include "split.hpp"
#include <SugarPP/include/sugarpp/range/enumerate.hpp>
#include <armadillo>
#include <boost/algorithm/string/predicate.hpp>
#include <map>

template<typename T>
void Phenotypes<T>::parse(std::istream &is) {
    int iid = 0;
    int phe = params.pheno_col;
    std::string line;
    arma::uword lineno = 0;
    std::set<std::string> sample_set;

    // column 0 is the original (*phenotypes)
    (*phenotypes).resize(params.nperms + 1);

    while (std::getline(is, line)) {
        if (boost::starts_with(line, "#") || lineno == 0) {// Skip the header
            lineno++;
            continue;
        }
        RJBUtil::Splitter<std::string_view> splitter(line, " \t");

        if (splitter[phe] == "NA" || splitter[phe] == "-9") {
            continue;
        } else {
            if (sample_set.find(splitter[iid]) != sample_set.end()) {
                throw(std::runtime_error(fmt::format("Duplicate sample ID: {}", splitter[iid])));
            }
            sample_set.insert(splitter[iid]);
            samples->push_back(splitter[iid]);
            (*phenotypes)[0].push_back(static_cast<int8_t>(std::stoi(splitter[phe])));
            lookup[splitter[iid]] = (*phenotypes)[0].back();
            // Checking for erroneous (*phenotypes)
            if ((*phenotypes)[0].back() < 0 || (*phenotypes)[0].back() > 1) {
                if (params.verbose) {
                    fmt::print(std::cerr, "{} {}\n", splitter[0], splitter[phe]);
                }
                throw(std::runtime_error(fmt::format("Incorrect phenotype value at lineno: {}", lineno)));
            }
            if (params.swap) {// Swap case-control status
                switch ((*phenotypes)[0].back()) {
                    case 1:
                        (*phenotypes)[0].back() = 0;
                        break;
                    case 0:
                        (*phenotypes)[0].back() = 1;
                        break;
                    default:
                        break;
                }
            }
        }
        lineno++;
    }
    IndexSort indexSort(*samples);
    indexSort.sort_vector(*samples);
    indexSort.sort_vector((*phenotypes)[0]);// Both must be sorted
    create_indexers();
    arma::wall_clock timer;
    timer.tic();
    // Read in the permutations from a TSV file where each row of the file is a permutation
    if (params.read_permutations) {
        std::ifstream ifs(*params.read_permutations);
        int i = 0;
        while (std::getline(ifs, line)) {
            RJBUtil::Splitter<std::string> splitter(line, " \t");
            if(splitter.size() == 0) {
                continue;
            }
            T perm;
            for (const auto &s : splitter) {
                perm.push_back(static_cast<int8_t>(std::stoi(s)));
            }
            (*phenotypes)[i + 1] = perm;
            i++;
        }
        ifs.close();
    } else {
        if (params.cov) {
            parse_cov();
            create_indexers();
        }
        shuffle();
    }
#if 0
    // Print the permutations to a file
    std::ofstream ofs("permutations.txt");
    for (const auto &p : *phenotypes) {
        for (const auto &[i, v] : Enumerate(p)) {
            if (i < p.size() - 1) {
                fmt::print(ofs, "{}\t", v);
            } else {
                fmt::print(ofs, "{}\n", v);
            }
        }
    }
    ofs.close();
#endif
    pad_phenotypes();
    fmt::print(std::cerr, "Time spent generating permutations: {}\n", timer.toc());
}

template<typename T>
void Phenotypes<T>::count_cov() {
    std::ifstream ifs(*params.cov);
    std::string line;
    unsigned long lineno = 0;
    while(std::getline(ifs, line)) {
        lineno++;
    }
    cov_lines = lineno - 1;
}

/**
 * @brief Parse covariate file
 * @tparam T Type of phenotype vector
 */
template<typename T>
void Phenotypes<T>::parse_cov() {
    using namespace RJBUtil;

    count_cov();

    std::ifstream ifs(*params.cov);
    std::string line;
    unsigned long lineno = 0;
    std::vector<std::string> cov_samples;

    // Parse covariate values
    std::map<unsigned long, std::vector<std::pair<unsigned long, std::string>>> failed_conversion;
    while(std::getline(ifs, line)) {
        lineno++;
        Splitter<std::string> split(line, " \t");
        if (lineno == 1) { // Handle header
            if (static_cast<long>(split.size()) - 1 <= 0) { // Error out for misformatted file
                throw(std::runtime_error("ERROR: Provided covariate file does not match expected format. sampleID y1 y2 y3 ..."));
            }
            cov = arma::mat(cov_lines, split.size() - 1, arma::fill::zeros);
            continue;
        }

        for (const auto &[col, v] : Enumerate(split)) {
            if (col == 0) {
                cov_samples.push_back(v);
                continue;
            }
            try {
                (*cov)(lineno - 2, col - 1) = std::stod(v);
            } catch (std::invalid_argument &e) {
                failed_conversion[col - 1].push_back({lineno, v});
            }
        }
    }
    // Convert failed conversions to a levels and then dummy variables
    std::vector<unsigned long> non_numeric_columns;
    std::map<unsigned long, std::unordered_map<std::string, unsigned long>> levels;
    for (const auto &[col, vec] : failed_conversion) {
        if (std::find(non_numeric_columns.begin(), non_numeric_columns.end(), col) == non_numeric_columns.end()) {
            non_numeric_columns.push_back(col);
        }
        for (const auto &[lineno, v] : vec) {
            levels[col].insert({v, 0});
        }
    }

    // Shed non-numeric columns
    auto nnc = arma::conv_to<arma::uvec>::from(non_numeric_columns);
    nnc.print("Non-numeric columns");
    for(int i = nnc.n_elem; i > 0; i--) {
        cov->shed_col(nnc(i - 1));
    }

    // Insert levels[col].size() - 1 dummy columns for each non-numeric column
    for (const auto &[col, level] : levels) {
        // Insert level.size() columns
        for (const auto &[i, l] : Enumerate(level)) {
            // remove the last column
            if(i == level.size() - 1) {
                levels[col].erase(l.first);
                break;
            }
            std::cerr << fmt::format("Adding column for level: {}\n", l.first);
            levels[col][l.first] = cov->n_cols;
            std::cerr << fmt::format("Column index: {}\n", l.second);
            cov->insert_cols(cov->n_cols, 1);
            std::cerr << fmt::format("Current number of columns: {}\n", cov->n_cols);
        }
        for (const auto &[lineno, v] : failed_conversion[col]) {
            if (levels[col].find(v) == levels[col].end()) {
                continue;
            } else {
                (*cov)(lineno - 2, levels[col][v]) = 1;
            }
        }
    }

    // Sort values
    IndexSort indexsort(cov_samples);
    indexsort.sort_vector(cov_samples);
    cov = (*cov)(arma::conv_to<arma::uvec>::from(indexsort.idx));

    // Verify that samples are present. If not, remove and count the number of cases and controls removed
    // Prune samples that lack covariates or phenotypes
    unsigned long cases_removed = 0;
    unsigned long controls_removed = 0;
    std::set<std::string> cov_set(cov_samples.begin(), cov_samples.end());
    std::set<std::string> pheno_set(samples->begin(), samples->end());
    if (cov_set.size() != pheno_set.size()) {
        std::vector<std::string> diff;
        std::set_symmetric_difference(cov_set.begin(), cov_set.end(), pheno_set.begin(), pheno_set.end(), std::back_inserter(diff));
        // Remove samples that are not in both files and count the number of cases and controls removed
        for (const auto &d : diff) {
            auto it = std::find(cov_samples.begin(), cov_samples.end(), d);
            if (it != cov_samples.end()) {
                auto idx = std::distance(cov_samples.begin(), it);
                cov_samples.erase(cov_samples.begin() + idx);
                (*cov).shed_row(idx);
                cases_removed += (*phenotypes)[0][idx];
                controls_removed += 1 - (*phenotypes)[0][idx];
            }
            it = std::find(samples->begin(), samples->end(), d);
            if (it != samples->end()) {
                auto idx = std::distance(samples->begin(), it);
                samples->erase(samples->begin() + idx);
                (*phenotypes)[0].erase((*phenotypes)[0].begin() + idx);
            }
        }
    }

    // Report the number of removed cases and controls
    fmt::print(std::cerr, "Cases removed: {}\n", cases_removed);
    fmt::print(std::cerr, "Controls removed: {}\n", controls_removed);

    // Check that all cov_samples equal samples
    for (const auto &[i, v] : Enumerate(cov_samples)) {
        if (v != (*samples)[i]) {
            throw(std::runtime_error("ERROR: Covariate file does not match phenotype file."));
        }
    }
}

/**
 * @brief Pad phenotypes to allow for vectorized reads past the end of the vector
 * @tparam T Type of phenotype vector
 *
 * Makes sure that we can do vectorized reads a little bit past the end of phenotypes without a segfault
 * Necessary for a vectorized gather from phenotypes, in Statistic::calculate
 */
template<typename T>
void Phenotypes<T>::pad_phenotypes() {
    for (auto &p : *phenotypes) {
        if (p.capacity() < p.size() + 3) {
            p.resize(p.size() + 3);
        }
    }
}

template<typename T>
void Phenotypes<T>::create_indexers() {
    int case_count = 0;
    int control_count = 0;
    int excluded = 0;
    for (const auto &p : phenotypes->at(0)) {
        switch (p) {
            case 1:
                case_count++;
                break;
            case 0:
                control_count++;
                break;
            default:
                excluded++;
                break;
        }
    }
    if (case_count <= 0 || control_count <= 0) {
        fmt::print(std::cerr, "ERROR: Case or Control count is equal to 0. Check your phenotype file.");
        std::exit(1);
    }

    fmt::print(std::cerr, "Phenotype counts --> cases: {}, controls: {}, excluded: {}\n", case_count, control_count, excluded);
    indexer = std::make_shared<Indexer<T>>(case_count, control_count, (*samples), (*phenotypes)[0]);
}

template<typename T>
void Phenotypes<T>::shuffle() {
    if (params.cov) {
        // Fit glm to covariates
        arma::mat design = arma::join_horiz(arma::ones<arma::mat>(cov->n_rows, 1), *cov);
        arma::vec Y((*phenotypes)[0].size());
        for (const auto &[i, v] : Enumerate((*phenotypes)[0])) {
            Y(i) = v;
        }
        Binomial link("logit");
        GLM glm(design, Y, link, params.algorithm);
        Permute<T> permute(params.seed);
        arma::vec odds = arma::exp(glm.eta_);

        permute.generate_permutations(phenotypes, odds, indexer->case_count, params.nperms, params.nthreads, params.epsilon);
    } else {
        for (int i = 1; i < params.nperms + 1; i++) {
            (*phenotypes)[i] = (*phenotypes)[i - 1];
            for (int j = (*phenotypes)[0].size() - 1; j > 0; j--) {
                std::uniform_int_distribution<> dis(0, j);
                int k = dis(gen);
                using std::swap;
                swap((*phenotypes)[i][j], (*phenotypes)[i][k]);
            }
        }
    }
}

template<typename T>
Phenotypes<T>::Phenotypes(Parameters params_, std::seed_seq &seed_source) : params(std::move(params_)), gen(seed_source) {
    samples = std::make_shared<std::vector<std::string>>();
    phenotypes = std::make_shared<std::vector<T>>();
    std::ifstream ifs(params.pheno);
    parse(ifs);
}

template class Phenotypes<pheno_vector>;
template class Phenotypes<compressed_pheno_vector>;
