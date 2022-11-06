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

    if (pairs.empty()) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            auto p = (*indexer).back_translate(it.row());
            pairs.emplace_back(p);
        }
    }

    for (auto &p : pairs) {
        auto &[p1, p2] = p;
        int8_t x = phenotypes_[p1];
        int8_t y = phenotypes_[p2];

        cscs += ((x == 1) && (y == 1));
        cscn += ((x == 1) && (y == 0));
        cscn += ((x == 0) && (y == 1));
        cncn += ((x == 0) && (y == 0));
    }

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