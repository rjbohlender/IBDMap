//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "math.hpp"

double cor(const pair_data_t &X, const pair_data_t &Y) {
    const double n = static_cast<double>(X.n_elem);
    double sum_x = 0;
    double sum_y = 0;
    double sum_x2 = 0;
    double sum_y2 = 0;
    double sum_xy = 0;

    for (auto it = X.begin(); it != X.end(); ++it) {
        const double value = static_cast<double>(*it);
        sum_x += value;
        sum_x2 += value * value;
    }
    for (auto it = Y.begin(); it != Y.end(); ++it) {
        const double value = static_cast<double>(*it);
        sum_y += value;
        sum_y2 += value * value;
    }

    auto x_it = X.begin();
    auto y_it = Y.begin();
    while (x_it != X.end() && y_it != Y.end()) {
        if (x_it.row() < y_it.row()) {
            ++x_it;
        } else if (y_it.row() < x_it.row()) {
            ++y_it;
        } else {
            sum_xy += static_cast<double>(*x_it) * static_cast<double>(*y_it);
            ++x_it;
            ++y_it;
        }
    }

    const double numerator = n * sum_xy - sum_x * sum_y;
    const double x_term = n * sum_x2 - sum_x * sum_x;
    const double y_term = n * sum_y2 - sum_y * sum_y;
    const double denominator = std::sqrt(x_term * y_term);
    const double r = numerator / denominator;
    return r * r;
};
