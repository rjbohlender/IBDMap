//
// Created by Bohlender,Ryan James on 2/11/21.
//

#include "inputvalidator.hpp"
#include "split.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <fmt/ostream.h>
#include <iostream>
#include <sys/stat.h>

bool exists(const std::string &path) {
    struct stat sb;
    return (stat(path.c_str(), &sb) == 0);
}

void InputValidator::check_gmap(const std::string &line, size_t line_no) {
    std::string line_trim = boost::trim_right_copy(line);
    RJBUtil::Splitter<std::string> splitter(line_trim, " \t");
    if (splitter.size() != 3) {
        fmt::print(std::cerr, "Incorrect number of columns at line {} in GMAP input.", line_no);
        std::exit(-1);
    }
    if (line_no == 0) {
        if (splitter[0] != "pos") {
            fmt::print(std::cerr, "Incorrect GMAP header. Columns should be pos chr cM, separated by tabs. The header is expected to have this format.");
            std::exit(-1);
        }
        if (splitter[1] != "chr") {
            fmt::print(std::cerr, "Incorrect GMAP header. Columns should be pos chr cM, separated by tabs. The header is expected to have this format.");
            std::exit(-1);
        }
        if (splitter[2] != "cM") {
            fmt::print(std::cerr, "Incorrect GMAP header. Columns should be pos chr cM, separated by tabs. The header is expected to have this format.");
            std::exit(-1);
        }
    } else {
        try {
            std::stoi(splitter[0]);
        } catch (std::exception &e) {
            fmt::print(std::cerr, "Failed to parse position at line {} in GMAP input.", line_no);
            std::exit(-1);
        }
        try {
            std::stod(splitter[2]);
        } catch (std::exception &e) {
            fmt::print(std::cerr, "Failed to parse CM at line {} in GMAP input.", line_no);
            std::exit(-1);
        }
    }
}

void InputValidator::check_pheno(const std::string &line, size_t line_no) {
    RJBUtil::Splitter<std::string> splitter(line, " \t");
    if (line_no == 0) {
        if (splitter.size() < 2) {
            fmt::print(std::cerr, "Phenotypes file must have at least 2 columns. The sampleid and the phenotype.");
            std::exit(-1);
        }
        pheno_column_count = splitter.size();
    }
    if (splitter.size() != pheno_column_count) {
        fmt::print(std::cerr, "Phenotypes file has an incorrect number of columns at line {}.", line_no);
        std::exit(-1);
    }
}
