//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#include "breakpoint.hpp"
#include "indexer.hpp"
#include "parameters.hpp"
#include "reporter.hpp"
#include "transposed_phenotypes.hpp"
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
    std::shared_ptr<TransposedPhenotypes> transposed;
    std::shared_ptr<BitPackedPhenotypes> bitpacked;
    uint64_t seq;
    Parameters params;
    std::shared_ptr<std::vector<T>> phenotypes;


    std::pair<std::vector<int32_t>, std::vector<int32_t>> pairs;
    std::vector<uint8_t> left_phenos;
    std::vector<uint8_t> right_phenos;

    int64_t lcm_common = 0;
    int64_t lcm_cscs_scale = 0;
    int64_t lcm_cscn_scale = 0;
    int64_t lcm_cncn_scale = 0;

    void setup_lcm();
    void initialize();
    void permute();
    void permute_bulk();
    void permute_bulk_bitpacked();
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
              uint64_t seq_,
              Parameters params_,
              std::shared_ptr<std::vector<T>> phenotypes_,
              std::shared_ptr<TransposedPhenotypes> transposed_ = nullptr,
              std::shared_ptr<BitPackedPhenotypes> bitpacked_ = nullptr);

    double calculate(T &phenotypes_, bool original_ = false);

    void run();
    void cleanup();
    static void test_statistic();

    friend void test_permutation_equivalence();
};

#endif//CARVAIBD_STATISTIC_HPP
