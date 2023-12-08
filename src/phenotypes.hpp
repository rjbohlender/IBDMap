//
// Created by Bohlender,Ryan J on 10/8/21.
//

#ifndef CARVAIBD_PHENOTYPES_HPP
#define CARVAIBD_PHENOTYPES_HPP

#include "indexer.hpp"
#include "parameters.hpp"
#include <istream>
#include <vector>
#include <memory>
#include <pcg_random.hpp>
#include <stocc/stocc.h>

template <typename T>
class Phenotypes {
    Parameters params;

    pcg64 gen;

    std::optional<arma::mat> cov;

    unsigned long cov_lines = 0;

    void parse(std::istream &is);
    void parse_cov();
    void count_cov();
    void create_indexers();
    void shuffle();
    void pad_phenotypes();
public:
    std::shared_ptr<std::vector<std::string>> samples;
    std::shared_ptr<std::vector<T>> phenotypes;
    std::shared_ptr<Indexer<T>> indexer;
    std::unordered_map<std::string, int8_t> lookup;

    explicit Phenotypes(Parameters params_, std::seed_seq &seed_source);

};


#endif//CARVAIBD_PHENOTYPES_HPP
