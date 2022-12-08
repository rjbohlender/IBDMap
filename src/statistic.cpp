//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include "statistic.hpp"
#include "split.hpp"
#include <fmt/include/fmt/ostream.h>
#include <utility>

#ifdef __AVX2__
#include <immintrin.h>
#endif

template <typename T>
Statistic<T>::Statistic(arma::SpCol<int32_t> data_,
                     Breakpoint bp_,
                     std::shared_ptr<Indexer<T>> indexer_,
                     std::shared_ptr<Reporter> reporter_,
                     Parameters params_,
                     std::shared_ptr<std::vector<T>> phenotypes_) : data(std::move(data_)), indexer(std::move(indexer_)),
                                                                               params(std::move(params_)),
                                                                               bp(std::move(bp_)), reporter(std::move(reporter_)),
                                                                               phenotypes(std::move(phenotypes_)) {
    if (params.enable_testing) {
        test_statistic();
    }
}

template <typename T>
double
Statistic<T>::calculate(T &phenotypes_, bool original_) noexcept {
    double statistic;

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    if (pairs.first.empty()) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            auto [left, right] = (*indexer).back_translate(it.row());
            pairs.first.emplace_back(left);
            pairs.second.emplace_back(right);
        }
    }

    size_t i = 0;

#if defined __AVX2__
    if constexpr (typeid(T) == typeid(int8_t)) {
        // AVX2 has less friendly instructions for sure. I wonder if I can clean up all these nasty casts.
        // We could avoid 2 expensive instructions if we stored phenotypes as 0xFF (or honestly just the most significant bit)
        // As is this, helps on Skylake-server but does not help on Desktop Zen 3, until the dataset gets too big for L2, then it helps.

        for (; i + 7 < pairs.first.size(); i += 8) {
            auto left_addresses = _mm256_loadu_si256((const __m256i *) &pairs.first[i]);
            auto right_addresses = _mm256_loadu_si256((const __m256i *) &pairs.second[i]);

            auto lefts = _mm256_i32gather_epi32((const int *) phenotypes_.data(), left_addresses, 1);
            auto rights = _mm256_i32gather_epi32((const int *) phenotypes_.data(), right_addresses, 1);

            const auto MASK_AND_SET_HIGH_BIT = _mm256_setr_epi8(0, 0, 0, 127, 0, 0, 0, 127, 0, 0, 0, 127, 0, 0, 0, 127,
                                                                0, 0, 0, 127, 0, 0, 0, 127, 0, 0, 0, 127, 0, 0, 0, 127);

            // Set the high bit on each byte only if the value was 1
            // We are using this both to mask out only the bytes we care about, as we retrieved 3 bytes of junk for every byte we want
            // And also to set things up for movemask_epi8 to look at only the most significant bit
            auto left_masked = _mm256_add_epi8(lefts, MASK_AND_SET_HIGH_BIT);
            auto right_masked = _mm256_add_epi8(rights, MASK_AND_SET_HIGH_BIT);

            auto left_packed = _mm256_movemask_epi8(left_masked);
            auto right_packed = _mm256_movemask_epi8(right_masked);

            auto cscs_batch = left_packed & right_packed;
            auto cscn_batch = left_packed ^ right_packed;

            cscs += __builtin_popcount(cscs_batch);
            cscn += __builtin_popcount(cscn_batch);
        }
    }
#endif

    for (; i < pairs.first.size(); ++i) {
        const auto x = phenotypes_[pairs.first[i]];
        const auto y = phenotypes_[pairs.second[i]];

        cscs += x & y;
        cscn += x ^ y;
    }

    cncn = pairs.first.size() - cscs - cscn;

    if (original_) {
        orig_cscs = static_cast<double>(cscs);
        orig_cscn = static_cast<double>(cscn);
        orig_cncn = static_cast<double>(cncn);
    }

    if (params.contcont) {
        statistic = static_cast<double>(cscs) / (*indexer).case_case - static_cast<double>(cscn) / (*indexer).case_cont - static_cast<double>(cncn) / (*indexer).cont_cont;
    } else {
        statistic = static_cast<double>(cscs) / (*indexer).case_case - static_cast<double>(cscn) / (*indexer).case_cont;
    }

    return statistic;
}

template <typename T>
void Statistic<T>::run() {
    initialize();

    arma::vec odds;

    permute();

    // Output
    std::stringstream ss;
    build_output(ss);
    if (params.verbose) {
        fmt::print(std::cerr, "Finished {}\t{}\n", bp.breakpoint.first, bp.breakpoint.second);
    }

    reporter->submit(ss.str());
    cleanup();
    done = true;
}

template <typename T>
void Statistic<T>::build_output(std::stringstream &ss) {
    double cscs = (*indexer).case_case;
    double cscn = (*indexer).case_cont;
    double cncn = (*indexer).cont_cont;
    fmt::print(ss, "{}\t{}\t", bp.breakpoint.first, bp.breakpoint.second);
    fmt::print(ss, "{}\t{}\t{}\t{}", orig_cscs / cscs, orig_cscn / cscn, orig_cncn / cncn, original);
    for (int i = 0; i < params.nperms; i++) {
        fmt::print(ss, "\t{}", permuted[i]);
    }
    fmt::print(ss, "\n");
}

template <typename T>
void Statistic<T>::permute() {
    double val;
    for (int i = 1; i <= params.nperms; i++) {
        val = calculate(phenotypes->at(i), false);
        permuted.push_back(val);
    }
}

template <typename T>
void Statistic<T>::initialize() {
    original = calculate(phenotypes->at(0), true);
}

template <typename T>
void Statistic<T>::cleanup() {
    data.reset();
    permuted.clear();
    permuted.shrink_to_fit();
}

template <typename T>
void Statistic<T>::test_statistic() {
    // Test case setup
    std::vector<std::string> tsamples{
            "case1", "case2", "case3", "case4", "case5",
            "case6", "case7", "case8", "case9", "case10",
            "control1", "control2", "control3", "control4", "control5",
            "control6", "control7", "control8", "control9", "control10"};
    T tphenotypes_{
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Indexer tindexer(10, 10, tsamples, tphenotypes_);

    arma::sp_vec tdata(tsamples.size() * (tsamples.size() - 1.) / 2.);

    // 4 case-case pairs, 4 control-control pairs, 4 case-control pairs
    tdata(tindexer.translate("case1", "case2")) = 1.;
    tdata(tindexer.translate("case1", "case3")) = 1.;
    tdata(tindexer.translate("case1", "case4")) = 1.;
    tdata(tindexer.translate("case5", "case6")) = 1.;
    tdata(tindexer.translate("control1", "control2")) = 1.;
    tdata(tindexer.translate("control1", "control3")) = 1.;
    tdata(tindexer.translate("control1", "control4")) = 1.;
    tdata(tindexer.translate("control5", "control6")) = 1.;
    tdata(tindexer.translate("case1", "control2")) = 1.;
    tdata(tindexer.translate("case1", "control3")) = 1.;
    tdata(tindexer.translate("case1", "control4")) = 1.;
    tdata(tindexer.translate("case5", "control6")) = 1.;

    // Follow normal execution
    double statistic;

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    std::vector<std::pair<size_t, size_t>> tpairs;

    if (tpairs.empty()) {
        for (auto it = tdata.begin(); it != tdata.end(); ++it) {
            auto p = tindexer.back_translate(it.row());
            try {
                tpairs.emplace_back(p);
            } catch (std::length_error &e) {
                std::cerr << "Failed to emplace or push at " << __LINE__ << std::endl;
                throw(e);
            }
        }
    }

    for (auto &p : tpairs) {
        auto &[p1, p2] = p;
        int x = tphenotypes_[p1];
        int y = tphenotypes_[p2];

        cscs += ((x == 1) && (y == 1));
        cscn += ((x == 1) && (y == 0));
        cscn += ((x == 0) && (y == 1));
        cncn += ((x == 0) && (y == 0));
    }

    statistic = static_cast<double>(cscs) / tindexer.case_case - static_cast<double>(cscn) / tindexer.case_cont - static_cast<double>(cncn) / tindexer.cont_cont;
    fmt::print(std::cerr, "Test statistic: {}\n", statistic);
    fmt::print(std::cerr, "cscs: {}, cscn: {}, cncn: {}\n", cscs, cscn, cncn);
    fmt::print(std::cerr,
               "cscs_count: {}, cscn_count: {}, cncn_count: {}\n",
               tindexer.case_case,
               tindexer.case_cont,
               tindexer.cont_cont);
}

template class Statistic<pheno_vector>;
template class Statistic<compressed_pheno_vector>;
