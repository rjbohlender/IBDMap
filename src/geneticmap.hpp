//
// Created by Bohlender,Ryan James on 10/22/19.
//

#ifndef CARVAIBD_GENETICMAP_HPP
#define CARVAIBD_GENETICMAP_HPP

#include "isgzipped.hpp"
#include "split.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <map>

class GeneticMap {
    std::map<std::string, std::map<int, double>> gmap_;

    /**
   * @brief Parser for the genetic map format.
   * @param is Input stream, possibly unzipped
   */
    void parse(std::istream &is) {
        // Format: pos	chr	cM -- tab separated, header in the file
        std::string line;
        while (std::getline(is, line)) {
            RJBUtil::Splitter<std::string_view> splitter(line, " \t");
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
        boost::iostreams::filtering_streambuf<boost::iostreams::input> streambuf;
        std::ifstream ifs;
        if (is_gzipped(path)) {
            ifs.open(path, std::ios_base::in | std::ios_base::binary);
            streambuf.push(boost::iostreams::gzip_decompressor());
            streambuf.push(ifs);
        } else {
            ifs.open(path, std::ios_base::in);
            streambuf.push(ifs);
        }
        std::istream is(&streambuf);
        parse(is);
    }

    /**
   * @brief Find the nearest bases for a given chromosome and position.
   * @param chrom Chromosome to search on
   * @param pos Position to search around
   * @return Pair of pairs with positions and recombination rates
   */
    std::pair<std::pair<int, double>, std::pair<int, double>> find_nearest(const std::string &chrom, int pos) {
        // If variant found in gmap then return pair with equal values, otherwise return prior and following.
        // The interpolation will drop out, leaving us with the correct value when the variant is found.
        auto lower = std::lower_bound(gmap_[chrom].begin(),
                                      gmap_[chrom].end(),
                                      pos,
                                      [](auto v1, int pos) { return v1.first < pos; });
        auto upper = std::upper_bound(gmap_[chrom].begin(),
                                      gmap_[chrom].end(),
                                      pos,
                                      [](int pos, auto v1) { return pos < v1.first; });
        if (lower == upper) {                   // Value not found
            if (lower == gmap_[chrom].begin()) {// Can't get lower value, return as is.
                return std::make_pair(*lower, *upper);
            } else if (lower == gmap_[chrom].end()) {// Position falls outside map, return highest value
                lower--;
                upper--;
                return std::make_pair(*lower, *upper);
            } else {// Get the value bracketing below
                lower--;
                return std::make_pair(*lower, *upper);
            }
        }
        // Value found
        upper--;// Point to the same position
        return std::make_pair(*lower, *upper);
    }
};

#endif//CARVAIBD_GENETICMAP_HPP
