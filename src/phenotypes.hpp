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

template <typename T>
class Phenotypes {
    Parameters params;

    std::mt19937_64 gen;


    void parse(std::istream &is);
    void create_indexers();
    void shuffle();
    void pad_phenotypes();
public:
    std::shared_ptr<std::vector<std::string>> samples;
    std::shared_ptr<std::vector<T>> phenotypes;
    std::shared_ptr<Indexer<T>> indexer;
    std::unordered_map<std::string, int8_t> lookup;

    explicit Phenotypes(Parameters params_);
};


#endif//CARVAIBD_PHENOTYPES_HPP
