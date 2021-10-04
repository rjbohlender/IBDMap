//
// Created by Bohlender,Ryan James on 9/24/19.
//

#ifndef CARVAIBD_PARAMETERS_HPP
#define CARVAIBD_PARAMETERS_HPP

#include <armadillo>
#include <fmt/include/fmt/ostream.h>
#include <optional>
#include <set>

/**
 * @brief Runtime parameters
 */
class Parameters {
public:
    std::string input;
    std::string pheno;
    std::string gmap;
    std::optional<std::string> info;
    size_t nperms;
    size_t nthreads;
    std::string output_path;
    unsigned seed = std::random_device{}();
    std::optional<arma::uword> lower_bound;
    bool swap = false;
    bool contcont = false;
    double min_dist;
    std::optional<double> rsquared;
    bool verbose = false;
    bool enable_testing = false;
    bool dash = false;
    std::optional<std::vector<int>> range;
    std::optional<std::vector<std::pair<int, int>>> exclude;
    std::optional<double> cM;
    std::optional<double> AF;
    std::optional<std::set<std::string>> sample_list;

    void print(std::ostream &os) {
        fmt::print(os, "Input: {}\n", input);
        fmt::print(os, "Pheno: {}\n", pheno);
        fmt::print(os, "Gmap: {}\n", gmap);
        if (info) {
            fmt::print(os, "Info: {}\n", *info);
        }
        fmt::print(os, "Output: {}\n", output_path);
        fmt::print(os, "Nperms: {}\n", nperms);
        fmt::print(os, "Nthreads: {}\n", nthreads);
        fmt::print(os, "Seed: {}\n", seed);
        if (lower_bound) {
            fmt::print(os, "Lower Bound: {}\n", *lower_bound);
        }
        fmt::print(os, "Swap: {}\n", swap);
        fmt::print(os, "Contcont: {}\n", contcont);
        fmt::print(os, "Min. Dist.: {}\n", min_dist);
        if (rsquared) {
            fmt::print(os, "R2: {}\n", *rsquared);
        }
        fmt::print(os, "Verbose: {}\n", verbose);
        fmt::print(os, "Enable Testing: {}\n", enable_testing);
        if (range) {
            fmt::print(os, "Range: {}\n", fmt::join(*range, "-"));
        }
        if (cM) {
            fmt::print(os, "cM: {}\n", *cM);
        }
        if (AF) {
            fmt::print(os, "AF: {}\n", *AF);
        }
        if (sample_list) {
            fmt::print(os, "Samples: {}\n", fmt::join(*sample_list, ","));
        }
    };
};
#endif//CARVAIBD_PARAMETERS_HPP
