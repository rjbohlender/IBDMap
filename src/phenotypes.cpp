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

        unsigned long col = 0;
        for (const auto &v : split) {
            if (col == 0) {
                cov_samples.push_back(v);
                col++;
                continue;
            }
            (*cov)(lineno - 2, col - 1) = std::stod(v);
            col++;
        }
    }

    // Sort values
    IndexSort indexsort(cov_samples);
    indexsort.sort_vector(cov_samples);
    cov = (*cov)(arma::conv_to<arma::uvec>::from(indexsort.idx));

    // Prune samples that lack covariates or phenotypes
    int i, j;
    if (cov_samples.size() > samples->size()) {
        for (i = 0, j = 0; i < cov_samples.size();) {
            if (cov_samples[i] == (*samples)[j]) {
                i++;
                j++;
            } else {
                while (cov_samples[i] != (*samples)[j]) {
                    cov_samples.erase(cov_samples.begin() + i);
                    (*cov).shed_row(i);
                }
            }
        }
    } else if (cov_samples.size() < samples->size()) {
        for (i = 0, j = 0; i < samples->size();) {
            if (cov_samples[i] == (*samples)[j]) {
                i++;
                j++;
            } else {
                while (cov_samples[i] != (*samples)[j]) {
                    samples->erase(samples->begin() + i);
                }
            }
        }
    }
    if (cov_samples.size() != samples->size()) {
        throw(std::runtime_error("ERROR: Covariate file does not match phenotype file."));
    }
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
    for (auto p : *phenotypes) {
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
        GLM glm(design, Y, link, "irls");
        Permute<T> permute(params.seed);
        arma::vec odds = arma::exp(glm.mu_);
        permute.generate_permutations(phenotypes, odds, indexer->case_count, params.nperms, params.nthreads, 0.0);
    } else {
        for (int i = 1; i < params.nperms + 1; i++) {
            (*phenotypes)[i] = (*phenotypes)[i - 1];
            for (int j = (*phenotypes)[0].size() - 1; j > 0; j--) {
                std::uniform_int_distribution<> dis(0, j);
                int k = dis(gen);
                std::swap((*phenotypes)[i][j], (*phenotypes)[i][k]);
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
