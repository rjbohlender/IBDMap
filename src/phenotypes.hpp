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

class Phenotypes {
    Parameters params;

    void parse(std::istream &is);
    void create_indexers();
    void build_groups(const std::map<std::vector<bool>, std::vector<arma::uword>> &fill_patterns);
public:
    std::shared_ptr<std::vector<std::string>> samples;
    std::vector<pheno_vector> phenotypes;
    std::shared_ptr<std::vector<Indexer>> indexer;
    std::optional<std::vector<std::vector<arma::uword>>> groups;

    explicit Phenotypes(Parameters params_);
};


#endif//CARVAIBD_PHENOTYPES_HPP
