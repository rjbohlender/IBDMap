//
// Created by Bohlender,Ryan James on 8/4/18.
//

#ifndef PERMUTE_ASSOCIATE_PERMUTATION_HPP
#define PERMUTE_ASSOCIATE_PERMUTATION_HPP

#include <armadillo>
#include <memory>
#include <stocc/stocc.h>
#include <thread>
#include <vector>

template<typename T>
struct Permute {
    Permute();
    explicit Permute(int seed);

    void generate_permutations(
            std::shared_ptr<std::vector<T>> permutations,
            arma::colvec &odds, arma::uword ncases, arma::uword nperm,
            arma::uword nthreads, double epsilon);

    void permute_thread(std::shared_ptr<std::vector<T>> p,
                        int ncases, int offset, int nperm, int seed);
    [[maybe_unused]] std::vector<T>
    epsilon_permutation(int nperm, arma::vec &odds, arma::uword ncases,
                        const std::string &transcript, double epsilon = 0.01);

    T unpack(int successes, int bin_size, bool shuffle,
             StochasticLib3 &rng);
    static void fisher_yates(T &x, StochasticLib3 &rng);

    StochasticLib3 sto;
    // Preserve group info for transcript
    bool bins_built;             // Whether the bins have been built or not -- for each
                                 // transcript
    std::vector<double> bin_odds;// Odds for the bins -- averaged within bin
    std::vector<int32_t> m;      // Number of samples in each bin
    std::vector<T> ret;
    arma::uvec sort_idx;// Indices of all samples in odds sorted order
    unsigned long long nsamples;
    void build_bins(arma::colvec &odds, double epsilon);
};
#endif// PERMUTE_ASSOCIATE_PERMUTATION_HPP
