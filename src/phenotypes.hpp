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

    std::mt19937_64 gen;

    void parse(std::istream &is);
    void create_indexers();
    void shuffle();
public:
    std::shared_ptr<std::vector<std::string>> samples;
    std::shared_ptr<std::vector<pheno_vector>> phenotypes;
    std::shared_ptr<Indexer> indexer;

    explicit Phenotypes(Parameters params_);
};


#endif//CARVAIBD_PHENOTYPES_HPP
