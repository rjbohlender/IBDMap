#include "statistic.hpp"
#include "math.hpp"

#include <cassert>
#include <cmath>

void test_permutation_equivalence();

int main() {
    pair_data_t x(5);
    pair_data_t y(5);
    x(1) = 1;
    x(2) = 2;
    x(4) = 4;
    y(0) = 1;
    y(2) = 2;
    y(3) = 3;

    const double expected_r2 = 121.0 / 476.0;
    assert(std::abs(cor(x, y) - expected_r2) < 1e-12);

    test_permutation_equivalence();
    return 0;
}
