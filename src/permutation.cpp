//
// Created by Bohlender,Ryan James on 8/4/18.
//

#include "permutation.hpp"
#include "jointhreads.hpp"
#include "types.hpp"
#include <cassert>

template<typename T>
Permute<T>::Permute() : sto(std::random_device{}()), bins_built(false) {}

template<typename T>
Permute<T>::Permute(int seed) : sto(seed), bins_built(false) {}

template<typename T>
void Permute<T>::generate_permutations(
        std::shared_ptr<std::vector<T>> permutations,
        arma::colvec &odds, arma::uword ncases, arma::uword nperm,
        arma::uword nthreads, double epsilon) {
    if (!bins_built) {
        build_bins(odds, epsilon);
    }
#ifndef NDEBUG
    arma::uword msize = std::accumulate(m.begin(), m.end(), 0);
    assert(msize == odds.n_rows);
#endif

    // Initialize permutations
    permutations->resize(nperm + 1);
    for (int i = 1; i < nperm + 1; i++) {
        (*permutations)[i].resize(odds.n_rows);
    }

    int step = nperm / nthreads;
    int remaining = nperm;
    {// Force scope so that threads automatically join and exit when done via
        // JoinThreads
        std::vector<std::thread> threads;
        JoinThreads thread_joiner(threads);
        for (int i = 0; i < nthreads; i++) {
            int seed = sto.IRandom(0, std::numeric_limits<int>::max());
            int offset = i * step + 1;
            if (remaining < 0) {
                std::cerr << "Failed during permutation.\n";
                std::exit(-1);
            }
            if (i == nthreads - 1) {
                threads.emplace_back(std::thread(&Permute::permute_thread, this,
                                                 permutations, ncases, offset,
                                                 remaining, seed));
            } else {
                threads.emplace_back(std::thread(&Permute::permute_thread, this,
                                                 permutations, ncases, offset, step,
                                                 seed));
            }
            remaining -= step;
        }
    }
}

template<typename T>
void Permute<T>::build_bins(arma::colvec &odds, double epsilon) {
    std::cerr << "Building bins." << std::endl;
    arma::wall_clock timer;
    timer.tic();
    sort_idx = arma::sort_index(odds);
    arma::vec odds_sorted = odds(sort_idx);

    double fudge = 1e-10;
    auto left = odds_sorted.begin();
    auto right = left;
    while (left != odds_sorted.end()) {
        double total = *left;
        while ((++right) != odds_sorted.end() && *right < (*left + epsilon + fudge)) {
            total += *right;
        }
        int dist = std::distance(left, right);
        double avg = total / dist;
        m.push_back(dist);
        bin_odds.push_back(avg);
        left = right;
    }

    std::cerr << "nbins: " << m.size() << " time: " << timer.toc() << std::endl;
    nsamples = sort_idx.n_elem;
    bins_built = true;
}

template<typename T>
void Permute<T>::permute_thread(
        std::shared_ptr<std::vector<T>> p, int ncases, int offset,
        int nperm, int seed) {
    StochasticLib3 rng(seed);
    std::vector<int32_t> tmp(bin_odds.size(), 0);

    for (int i = 0; i < nperm; i++) {
        rng.MultiFishersNCHyp(&tmp[0], &(m[0]), &(bin_odds[0]), ncases,
                              bin_odds.size());
        // Unpack bins
        arma::uword filled = 0;
        for (int j = 0; j < m.size(); j++) {// for each bin
            T r = unpack(
                    tmp[j], m[j], true, rng);// unpack and randomize cases into a vector
            for (int k = 0; k < r.size(); k++) {
                (*p)[offset + i][sort_idx[k + filled]] = r[k];
            }
            filled += m[j];
        }
    }
}

/**
 * @brief Unpack the permuted values into random order for a bin
 * @param successes Number of cases within bin
 * @param bin_size Total number of samples in bin
 * @param shuffle Whether to randomize or not
 * @return Vector of phenotype states
 */
template<typename T>
T Permute<T>::unpack(int successes, int bin_size, bool shuffle,
                     StochasticLib3 &rng) {
    T r(bin_size, 0);
    if (successes > 0) {
        for (int i = 0; i < successes; i++) {
            r[i] = 1;
        }
    }
    // Fisher-Yates Shuffle
    if (shuffle) {
        for (int i = r.size() - 1; i >= 1; --i) {
            auto j = rng.IRandom(0, i);
            using std::swap;
            swap(r[i], r[j]);
        }
    }
    assert(std::accumulate(r.begin(), r.end(), 0) == successes);
    return r;
}

template<typename T>
[[maybe_unused]] std::vector<T>
Permute<T>::epsilon_permutation(int nperm, arma::vec &odds, arma::uword ncases,
                                const std::string &transcript, double epsilon) {
    if (!bins_built) {
        build_bins(odds, epsilon);
    }
#ifndef NDEBUG
    arma::uword msize = std::accumulate(m.begin(), m.end(), 0);
    assert(msize == odds.n_rows);
#endif
    for (int i = 0; i < nperm; i++) {
        std::vector<int32_t> tmp(bin_odds.size(), 0);
        sto.MultiFishersNCHyp(&tmp[0], &(m[0]), &(bin_odds[0]), ncases,
                              bin_odds.size());

        // Unpack bins
        arma::uword filled = 0;
        for (int j = 0; j < m.size(); j++) {// for each bin
            T r = unpack(
                    tmp[j], m[j], true, sto);   // unpack and randomize cases into a vector
            for (int k = 0; k < r.size(); k++) {//
                ret[i][sort_idx[k + filled]] = r[k];
            }
            filled += m[j];
        }
    }
    return ret;
}

template<typename T>
void Permute<T>::fisher_yates(T &x, StochasticLib3 &rng) {
    for (int i = x.size() - 1; i >= 1; --i) {
        auto j = rng.IRandom(0, i);
        arma::uword tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
    }
}
template class Permute<pheno_vector>;
template class Permute<compressed_pheno_vector>;
