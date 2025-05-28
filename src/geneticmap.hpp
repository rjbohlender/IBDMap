//
// Created by Bohlender,Ryan James on 10/22/19.
//

#ifndef CARVAIBD_GENETICMAP_HPP
#define CARVAIBD_GENETICMAP_HPP

#include "inputvalidator.hpp"
#include "split.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/python.hpp>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <zstr.hpp>

using namespace RJBUtil;
namespace p = boost::python;

class GeneticMap {
    std::map<std::string, std::map<int, double>> gmap_;

    /**
   * @brief Parser for the genetic map format.
   * @param is Input stream, possibly unzipped
   */
    void parse(std::istream &is) {
        // Format: pos	chr	cM -- tab separated, header in the file
        std::string line;
        int lineno = -1;
        while (std::getline(is, line)) {
            lineno++;
            InputValidator::check_gmap(line, lineno);
            Splitter<std::string_view> splitter(line, " \t");
            if (splitter[0] == "pos") {// Skip header
                continue;
            }
            std::string chrom;
            int pos;
            double gpos;
            if (!boost::starts_with(splitter[1], "chr")) {
                std::stringstream ss;
                ss << "chr" << splitter[1];
                chrom = ss.str();
            } else {
                chrom = splitter[1];
            }
            pos = std::stoi(splitter[0]);
            gpos = std::stod(splitter[2]);
            if (gmap_.find(chrom) == gmap_.end()) {
                gmap_[chrom] = {};
                gmap_[chrom][pos] = gpos;
            } else {
                gmap_[chrom][pos] = gpos;
            }
        }
    }

public:
    /**
   * @brief Collects multiple genetic maps into a single object
   * @param paths Set of paths, one or more
   */
    explicit GeneticMap(const std::string &path) {
        if (path.empty()) {
            throw std::runtime_error("ERROR: no filepath given to the genetic map object.");
        }
        zstr::ifstream ifs(path);
        parse(ifs);
    }

    explicit GeneticMap(p::list paths) {
        if (p::len(paths) == 0) {
            throw std::runtime_error("ERROR: no filepath given to the genetic map object.");
        }
        while(p::len(paths) > 0) {
            std::string p = p::extract<std::string>(paths.pop());
            zstr::ifstream ifs(p);
            parse(ifs);
        }
    }

    /**
   * @brief Find the nearest bases for a given chromosome and position.
   * @param chrom Chromosome to search on
   * @param pos Position to search around
   * @return Pair of pairs with positions and recombination rates
   */
    [[nodiscard]] std::pair<std::pair<int, double>, std::pair<int, double>> find_nearest(const std::string &chrom, int pos) const {
        // If variant found in gmap then return pair with equal values, otherwise return prior and following.
        // The interpolation will drop out, leaving us with the correct value when the variant is found.
        auto lower = std::lower_bound(gmap_.at(chrom).begin(),
                                      gmap_.at(chrom).end(),
                                      pos,
                                      [](auto v1, int pos) { return v1.first < pos; });
        if (lower->first != pos) {                   // Value not found
            if (lower == gmap_.at(chrom).begin()) {// Can't get lower value, return as is.
                auto upper = lower;
                upper++;
                return std::make_pair(*lower, *upper);
            } else if (lower == gmap_.at(chrom).end()) {// Position falls outside map, return highest value
                lower--;
                auto upper = lower;
                return std::make_pair(*lower, *upper);
            } else {// Get the value bracketing below
                auto upper = lower;
                lower--;
                return std::make_pair(*lower, *upper);
            }
        }
        // Value found
        auto upper = lower;;// Point to the same position
        return std::make_pair(*lower, *upper);
    }

    /**
     * @brief Check if a given position on a chromosome is in the dictionary.
     * @param chrom The chromosome id
     * @param pos The position in basepairs
     * @return bool indicating presence or absence.
     */
    [[nodiscard]] bool contains(const std::string &chrom, int pos) const {
      return gmap_.at(chrom).count(pos) > 0;
    }
};

#endif//CARVAIBD_GENETICMAP_HPP
