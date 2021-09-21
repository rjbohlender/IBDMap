//
// Created by Bohlender,Ryan James on 9/22/20.
//

#ifndef CARVAIBD_INFO_HPP
#define CARVAIBD_INFO_HPP

#include "parameters.hpp"
#include <map>
#include <string>
#include <vector>

class Info {
    std::map<std::string, std::vector<double>> data;
    std::map<std::string, int> field_map;

public:
    explicit Info(std::istream &ifs);

    double get_field(const std::string &segment, const std::string &field);
    bool filter_segment(const std::string &segment, const Parameters &params);
};

#endif//CARVAIBD_INFO_HPP
