#include <cstdio>
#include <iostream>
#include <random>

#include <benchmark/benchmark.h>

static void BM_StringCreation(benchmark::State &state) {
    for (auto _ : state)
        std::string empty_string;
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation);

// Define another benchmark
static void BM_StringCopy(benchmark::State &state) {
    std::string x = "hello";
    for (auto _ : state)
        std::string copy(x);
}
BENCHMARK(BM_StringCopy);

static void BM_fp(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs{
            {0, 1},
            {0, 2},
            {0, 3},
            {0, 4},
            {0,6},
            {0,8},
            {1, 3},
            {1,5},
            {1,7},
            {1,9},
            {2,3},
            {2,5},
            {2,6},
            {2,8},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
    };
    std::vector<int8_t> phenotypes_{0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1};

    double cscs = 0;
    double cscn = 0;
    double cncn = 0;

    for (auto _ : state) {
        for (auto &p : pairs) {
            auto &[p1, p2] = p;
            int x = phenotypes_[p1];
            int y = phenotypes_[p2];

            if (x < -1 || x > 1 || y < -1 || y > 1) {
                throw(std::runtime_error("ERROR: invalid phenotype in calculate."));
            }

            cscs += ((x == 1) && (y == 1));
            cscn += ((x == 1) && (y == 0));
            cscn += ((x == 0) && (y == 1));
            cncn += ((x == 0) && (y == 0));
        }
    }

    benchmark::DoNotOptimize(cscs + cscn + cncn);
}

BENCHMARK(BM_fp);

static void BM_nochecks_fp(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs{
            {0, 1},
            {0, 2},
            {0, 3},
            {0, 4},
            {0,6},
            {0,8},
            {1, 3},
            {1,5},
            {1,7},
            {1,9},
            {2,3},
            {2,5},
            {2,6},
            {2,8},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
    };
    std::vector<int8_t> phenotypes_{0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1};

    double cscs = 0;
    double cscn = 0;
    double cncn = 0;

    for (auto _ : state) {
        for (auto &p : pairs) {
            auto &[p1, p2] = p;
            int x = phenotypes_[p1];
            int y = phenotypes_[p2];

            cscs += ((x == 1) && (y == 1));
            cscn += ((x == 1) && (y == 0));
            cscn += ((x == 0) && (y == 1));
            cncn += ((x == 0) && (y == 0));
        }
    }

    benchmark::DoNotOptimize(cscs + cscn + cncn);
}

BENCHMARK(BM_nochecks_fp);

static void BM_int(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs{
            {0, 1},
            {0, 2},
            {0, 3},
            {0, 4},
            {0,6},
            {0,8},
            {1, 3},
            {1,5},
            {1,7},
            {1,9},
            {2,3},
            {2,5},
            {2,6},
            {2,8},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
    };
    std::vector<int8_t> phenotypes_{0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1};

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    for (auto _ : state) {
        for (auto &p : pairs) {
            auto &[p1, p2] = p;
            auto x = phenotypes_[p1];
            auto y = phenotypes_[p2];

            if (x < -1 || x > 1 || y < -1 || y > 1) {
                throw(std::runtime_error("ERROR: invalid phenotype in calculate."));
            }

            cscs += ((x == 1) && (y == 1));
            cscn += ((x == 1) && (y == 0));
            cscn += ((x == 0) && (y == 1));
            cncn += ((x == 0) && (y == 0));
        }
    }

    benchmark::DoNotOptimize(cscs + cscn + cncn);
}

BENCHMARK(BM_int);

static void BM_nochecks_int(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs{
            {0, 1},
            {0, 2},
            {0, 3},
            {0, 4},
            {0,6},
            {0,8},
            {1, 3},
            {1,5},
            {1,7},
            {1,9},
            {2,3},
            {2,5},
            {2,6},
            {2,8},
            {3, 4},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 8},
    };
    std::vector<int8_t> phenotypes_{0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1};

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    for (auto _ : state) {
        for (auto &p : pairs) {
            auto &[p1, p2] = p;
            auto x = phenotypes_[p1];
            auto y = phenotypes_[p2];

            cscs += ((x == 1) && (y == 1));
            cscn += ((x == 1) && (y == 0));
            cscn += ((x == 0) && (y == 1));
            cncn += ((x == 0) && (y == 0));
        }
    }

    benchmark::DoNotOptimize(cscs + cscn + cncn);
}

BENCHMARK(BM_nochecks_int);

BENCHMARK_MAIN();
