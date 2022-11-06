//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#define ARMA_DONT_USE_WRAPPER

#include "breakpoint.hpp"
#include "indexer.hpp"
#include "parameters.hpp"
#include "reporter.hpp"
#include "types.hpp"
#include <armadillo>
#include <boost/optional.hpp>
#include <future>
#include <string>

/**
 * @brief IBD Statistic Calculation Class
 */
class Statistic {
    arma::sp_colvec data;
    std::shared_ptr<Indexer> indexer;
    std::shared_ptr<Reporter> reporter;
    Parameters params;
    std::shared_ptr<std::vector<pheno_vector>> phenotypes;


    std::vector<std::pair<arma::sword, arma::sword>> pairs;

    void initialize();
    void permute();
    void build_output(std::stringstream &ss);

public:
    bool done = false;

    Breakpoint bp;

    double original;
    double orig_cscs;
    double orig_cscn;
    double orig_cncn;

    std::vector<double> permuted;

    Statistic(arma::sp_colvec data_,
              Breakpoint bp_,
              std::shared_ptr<Indexer> indexer_,
              std::shared_ptr<Reporter> reporter_,
              Parameters params_,
              std::shared_ptr<std::vector<pheno_vector>> phenotypes_);

    double calculate(pheno_vector &phenotypes_, bool original_ = false) noexcept;

    void run();
    void cleanup();
    static void test_statistic();
};

#endif//CARVAIBD_STATISTIC_HPP
