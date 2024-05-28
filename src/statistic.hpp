//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

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
template <typename T>
class Statistic {
    arma::SpCol<int32_t> data;
    std::shared_ptr<Indexer<T>> indexer;
    std::shared_ptr<Reporter> reporter;
    Parameters params;
    std::shared_ptr<std::vector<T>> phenotypes;


    std::pair<std::vector<int32_t>, std::vector<int32_t>> pairs;

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

    Statistic(arma::SpCol<int32_t> data_,
              Breakpoint bp_,
              std::shared_ptr<Indexer<T>> indexer_,
              std::shared_ptr<Reporter> reporter_,
              Parameters params_,
              std::shared_ptr<std::vector<T>> phenotypes_);

    double calculate(T &phenotypes_, bool original_ = false) noexcept;

    void run();
    void cleanup();
    static void test_statistic();
};

#endif//CARVAIBD_STATISTIC_HPP
