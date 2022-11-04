//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#define ARMA_DONT_USE_WRAPPER

#include "types.hpp"
#include "breakpoint.hpp"
#include "indexer.hpp"
#include "parameters.hpp"
#include "reporter.hpp"
#include <string>
#include <armadillo>
#include <boost/optional.hpp>
#include <future>

/**
 * @brief IBD Statistic Calculation Class
 */
class Statistic {
    arma::sp_colvec data;
    std::shared_ptr<std::vector<Indexer>> indexer;
    std::vector<pheno_vector> phenotypes;
    Parameters params;

    std::shared_ptr<Reporter> reporter;

    std::vector<std::pair<arma::sword, arma::sword>> pairs;

    std::mt19937_64 gen;

    static void x1(int y, double &cscs, double &cscn);
    static void x0(int y, double &cscn, double &cncn);
    void initialize();
    static void joint_shuffle(std::vector<pheno_vector> &phen, std::mt19937_64 &gen);
    void joint_permute();
    static void group_unpack(std::vector<pheno_vector> &p_original,
                             const std::vector<pheno_vector> &p_tmp,
                             const std::vector<arma::uword> &group_indices);
    static void group_pack(const std::vector<pheno_vector> &p_original,
                           std::vector<pheno_vector> &p_tmp,
                           const std::vector<arma::uword> &groupIndices);
    void build_output(std::stringstream &ss);

public:
    bool done = false;

    Breakpoint bp;

    std::optional<std::vector<std::vector<arma::uword>>> groups;

    std::vector<double> original;
    std::vector<double> successes;
    std::vector<size_t> permutations;

    std::vector<std::vector<double>> permuted;
    std::vector<std::vector<double>> permuted_cscs;
    std::vector<std::vector<double>> permuted_cscn;
    std::vector<std::vector<double>> permuted_cncn;

    Statistic(arma::sp_colvec data_, Breakpoint bp_, std::shared_ptr<std::vector<Indexer>> indexer_, std::shared_ptr<Reporter> reporter_, Parameters params_, std::optional<std::vector<std::vector<arma::uword>>> groups_);

    double calculate(pheno_vector &phenotypes_, double cscs_count, double cscn_count, double cncn_count, size_t k);

    void run();
    void cleanup();
    static void test_statistic();
    static void test_group_permutation(unsigned int seed);
};

#endif//CARVAIBD_STATISTIC_HPP
