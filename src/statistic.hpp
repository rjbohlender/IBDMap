//
// Created by Bohlender,Ryan James on 9/23/19.
//

#ifndef CARVAIBD_STATISTIC_HPP
#define CARVAIBD_STATISTIC_HPP

#include <armadillo>
#include <future>
#include "parser.hpp"
#include "reporter.hpp"
#include "parameters.hpp"

/**
 * @brief IBD Statistic Calculation Class
 */
class Statistic {
    arma::sp_colvec data;
    std::vector<Indexer> &indexer;
    Parser<std::string> &parser;
    Parameters &params;

    std::shared_ptr<Reporter> reporter;

    std::vector<std::pair<arma::sword, arma::sword>> pairs;
public:
    bool done = false;

    Breakpoint &bp;

    std::vector<double> original;
    std::vector<double> successes;
    std::vector<double> permutations;

    std::vector<std::vector<double>> permuted;
    std::vector<double> permuted_cscs;
    std::vector<double> permuted_cscn;

    std::vector<arma::sword> rows;

    Statistic(arma::sp_colvec &&data_,
              Breakpoint &bp_,
              std::vector<Indexer> &indexer_,
              Parser<std::string> &parser_,
              std::shared_ptr<Reporter> reporter_,
              Parameters &params_);

    double calculate(std::vector<int> &phenotypes_, double case_case, double case_cont, double cont_cont);

    void run();
};

#endif //CARVAIBD_STATISTIC_HPP
