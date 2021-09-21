//
// Created by Bohlender,Ryan James on 12/11/19.
//

#ifndef CARVAIBD_INDEXSORT_HPP
#define CARVAIBD_INDEXSORT_HPP

#include <algorithm>// std::sort
#include <iostream>
#include <numeric>// std::iota
#include <vector>

/**
 * @brief Sort an index vector by the values of an associated data vector
 */
template<typename T>
class IndexSort {
public:
    std::vector<std::size_t> idx;

    explicit IndexSort(const std::vector<T> &v) {
        // initialize original index locations
        idx = std::vector<std::size_t>(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        sort(idx.begin(), idx.end(),
             [&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });
    }

    template<typename X>
    void sort_vector(std::vector<X> &v) {
        if (v.size() != idx.size()) {
            throw(std::runtime_error("Trying to use indices to sort a vector of different size."));
        }
        std::vector<X> r = v;
        for (std::size_t i = 0; i < v.size(); i++) {// Swap elements
            if (idx[i] == i) {
                continue;
            }
            r[i] = v[idx[i]];
        }
        v = r;
    }
};

#endif//CARVAIBD_INDEXSORT_HPP
