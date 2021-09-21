//
// Created by Bohlender,Ryan James on 12/13/19.
//

#ifndef CARVAIBD_SRC_BREAKPOINT_HPP
#define CARVAIBD_SRC_BREAKPOINT_HPP

#include <utility>
#include <vector>

/**
 * @brief Basic POD type to hold breakpoint information
 */
struct Breakpoint {
    std::pair<std::string, std::string> breakpoint;
    std::vector<double> segment_lengths;
    std::vector<std::pair<std::string, std::string>> ibd_pairs;
};

#endif//CARVAIBD_SRC_BREAKPOINT_HPP
