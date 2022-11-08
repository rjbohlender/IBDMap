//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include "statistic.hpp"
#include "split.hpp"
#include <SugarPP/include/sugarpp/range/enumerate.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <fmt/include/fmt/ostream.h>
#include <string>
#include <utility>
#include <valarray>

#ifdef __AVX512F__
#include <immintrin.h>

// Maybe I could use c++20's std::popcount instead??
static inline int32_t popcnt128(__m128i n) {
    const __m128i n_hi = _mm_unpackhi_epi64(n, n);
    return __builtin_popcountll(_mm_cvtsi128_si64(n)) + __builtin_popcountll(_mm_cvtsi128_si64(n_hi));
}
#endif

Statistic::Statistic(arma::sp_colvec data_,
                     Breakpoint bp_,
                     std::shared_ptr<Indexer> indexer_,
                     std::shared_ptr<Reporter> reporter_,
                     Parameters params_,
                     std::shared_ptr<std::vector<pheno_vector>> phenotypes_) : data(std::move(data_)), indexer(std::move(indexer_)),
                                                                               params(std::move(params_)),
                                                                               bp(std::move(bp_)), reporter(std::move(reporter_)),
                                                                               phenotypes(std::move(phenotypes_)) {
    if (params.enable_testing) {
        test_statistic();
    }
}

double
Statistic::calculate(pheno_vector &phenotypes_, bool original_) noexcept {
    double statistic;

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    int64_t cscs2 = 0;
    int64_t cscn2 = 0;

    if (pairs.first.empty()) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            auto [left, right] = (*indexer).back_translate(it.row());
            try {
                pairs.first.emplace_back(left);
                pairs.second.emplace_back(right);
            } catch (std::length_error &e) {
                fmt::print(std::cerr, "Failed to emplace or push at line {}.", __LINE__);
                throw(e);
            }
        }
    }

    size_t i = 0;

    /**
#ifdef __AVX512F__
    // As much as I want to use this, it's much faster for this specific method, but slows the whole application by a hair
    // On Ice Lake or newer, this is probably the winner!
    if  (phenotypes_.capacity() < phenotypes_.size() + 3) {
        phenotypes_.resize(phenotypes_.size() + 3);
    }

    for (; i + 15 < pairs.first.size(); i+= 16) {
        auto left_addresses = _mm512_loadu_si512(&pairs.first[i]);
        auto right_addresses = _mm512_loadu_si512(&pairs.second[i]);

        auto lefts = _mm512_i32gather_epi32(left_addresses, phenotypes_.data(), 1);
        auto rights = _mm512_i32gather_epi32(right_addresses, phenotypes_.data(), 1);
        
        auto left_packed = _mm512_cvtepi32_epi8(lefts);
        auto right_packed = _mm512_cvtepi32_epi8(rights);

        auto cscs_batch = _mm_and_si128(left_packed, right_packed);
        auto cscn_batch = _mm_xor_si128(left_packed, right_packed);

        cscs += popcnt128(cscs_batch);
        cscn += popcnt128(cscn_batch);
    }
#endif
     */

    for (; i + 1 < pairs.first.size(); i += 2) {

        const auto x1 = phenotypes_[pairs.first[i]];
        const auto y1 = phenotypes_[pairs.second[i]];

        cscs += x1 & y1;
        cscn += x1 ^ y1;

        auto x2 = phenotypes_[pairs.first[i + 1]];
        auto y2 = phenotypes_[pairs.second[i + 1]];

        cscs2 += x2 & y2;
        cscn2 += x2 ^ y2;
    }

    for (; i < pairs.first.size(); ++i) {
        const auto x = phenotypes_[pairs.first[i]];
        const auto y = phenotypes_[pairs.second[i]];

        cscs += x & y;
        cscn += x ^ y;
    }

    cscs += cscs2;
    cscn += cscn2;

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

void Statistic::run() {
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

void Statistic::build_output(std::stringstream &ss) {
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

void Statistic::permute() {
    double val;
    for (int i = 1; i <= params.nperms; i++) {
        val = calculate(phenotypes->at(i), false);
        permuted.push_back(val);
    }
}

void Statistic::initialize() {
    permuted.emplace_back();
    original = calculate(phenotypes->at(0), true);
}

void Statistic::cleanup() {
    data.reset();
    permuted.clear();
    permuted.shrink_to_fit();
}

void Statistic::test_statistic() {
    // Test case setup
    std::vector<std::string> tsamples{
            "case1", "case2", "case3", "case4", "case5",
            "case6", "case7", "case8", "case9", "case10",
            "control1", "control2", "control3", "control4", "control5",
            "control6", "control7", "control8", "control9", "control10"};
    pheno_vector tphenotypes_{
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