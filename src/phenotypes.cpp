//
// Created by Bohlender,Ryan J on 10/8/21.
//

#include "phenotypes.hpp"
#include "split.hpp"
#include <SugarPP/include/sugarpp/range/enumerate.hpp>
#include <armadillo>
#include <boost/algorithm/string/predicate.hpp>
#include <map>

void Phenotypes::parse(std::istream &is) {
    int iid = 0;
    int phe = 1;
    std::string line;
    arma::uword lineno = 0;
    std::map<std::vector<bool>, std::vector<arma::uword>> fill_patterns;

    while (std::getline(is, line)) {
        if (boost::starts_with(line, "#") || lineno == 0) {// Skip the header
            lineno++;
            continue;
        }
        RJBUtil::Splitter<std::string_view> splitter(line, " \t");

        samples->push_back(splitter[iid]);
        std::vector<bool> pattern;
        for (int i = phe; i < splitter.size(); i++) {
            if (phenotypes.size() < i) {
                phenotypes.emplace_back();
            }
            if (splitter[i] == "NA") {
                pattern.push_back(false);
                phenotypes[i - 1].push_back(-1);
            } else {
                pattern.push_back(true);
                phenotypes[i - 1].push_back(std::stoi(splitter[i]));
                // Checking for erroneous phenotypes
                if (phenotypes[i - 1].back() < 0 || phenotypes[i - 1].back() > 1) {
                    if (params.verbose) {
                        fmt::print(std::cerr, "{} {}\n", splitter[0], splitter[i]);
                    }
                    throw(std::runtime_error(fmt::format("Incorrect phenotype value at lineno: {}", lineno)));
                }
                if (params.swap) {// Swap case-control status
                    switch (phenotypes[i - 1].back()) {
                        case 1:
                            phenotypes[i - 1].back() = 0;
                            break;
                        case 0:
                            phenotypes[i - 1].back() = 1;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
        if (fill_patterns.count(pattern) == 0) {
            fill_patterns.emplace(std::make_pair(pattern, std::vector<arma::uword>({lineno - 1})));
        } else {
            fill_patterns[pattern].push_back(lineno - 1);
        }
        lineno++;
    }
    create_indexers();
    build_groups(fill_patterns);
}

void Phenotypes::create_indexers() {
    for (auto [i, p] : Enumerate(phenotypes)) {
        int case_count = 0;
        int control_count = 0;
        int excluded = 0;
        for (const auto &v : p) {
            switch (v) {
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
        if (params.verbose) {
            fmt::print(std::cerr, "Phenotype {}, cases: {}, controls: {}, excluded: {}\n", i, case_count, control_count, excluded);
        }
        indexer->emplace_back(Indexer(case_count, control_count, (*samples), p));
    }
}

void Phenotypes::build_groups(const std::map<std::vector<bool>, std::vector<arma::uword>> &fill_patterns) {
    if (fill_patterns.size() > 1) {
        groups = std::vector<std::vector<arma::uword>>();
        for (const auto &v : fill_patterns) {
            groups->push_back(v.second);
        }
        if (params.verbose) {
            std::cerr << "Groups: " << groups->size() << std::endl;
            std::cerr << "Group sizes: ";
        }
        for (const auto &v : fill_patterns) {
            if (params.verbose) {
                std::cerr << v.second.size() << " ";
            }
        }
        if (params.verbose) {
            std::cerr << std::endl;
        }
    }
}

Phenotypes::Phenotypes(Parameters params_) : params(std::move(params_)) {
    samples = std::make_shared<std::vector<std::string>>();
    indexer = std::make_shared<std::vector<Indexer>>();
    std::ifstream ifs(params.pheno);
    parse(ifs);
}
