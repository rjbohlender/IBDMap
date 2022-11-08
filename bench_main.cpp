#include <cstdio>
#include <iostream>
#include <random>
#include <valarray>
#include <array>

#include <immintrin.h>

#include <benchmark/benchmark.h>

static const auto NUM_PHENOS = 500000;
static const auto NUM_PAIRS = 256;

static std::vector<int8_t> make_random_bools(size_t n = NUM_PHENOS)
{
    auto seed = std::random_device{}();
    std::mt19937 rng(seed);
    std::uint32_t data;
    std::vector<int8_t> bools;
    bools.reserve(n);

    int bit_count = 0;

    for (size_t i = 0 ; i < n ; ++i)
    {
        if (bit_count == 0)
        {
            data = rng();
            bit_count = 32;
        }
        bool bit = data & 1;
        data >>= 1;
        bit_count--;
        bools.emplace_back(bit);
    }

    return bools;
}

template<typename T = size_t>
static std::vector<std::pair<T, T>> make_pairs(size_t n_phenos = NUM_PHENOS, size_t n_pairs = NUM_PAIRS) {
    auto seed = std::random_device{}();
    std::mt19937 rng(seed);
    std::uniform_int_distribution<T> dis(0, n_phenos - 1);
    std::vector<std::pair<T, T>> pairs;
    pairs.reserve(n_pairs);

    for (size_t i = 0; i < n_pairs; i++) {
        pairs.emplace_back(std::make_pair(dis(rng), dis(rng)));
    }

    std::sort(pairs.begin(), pairs.end());

    return pairs;
}

template<typename T = size_t>
static std::pair<std::vector<T>, std::vector<T>> make_pairs_struct_of_arrays(size_t n_phenos = NUM_PHENOS, size_t n_pairs = NUM_PAIRS) {
    auto seed = std::random_device{}();
    std::mt19937 rng(seed);
    std::uniform_int_distribution<T> dis(0, n_phenos - 1);
    std::vector<std::pair<T, T>> pairs;
    pairs.reserve(n_pairs);

    for (size_t i = 0; i < n_pairs; i++) {
        pairs.emplace_back(std::make_pair(dis(rng), dis(rng)));
    }

    std::sort(pairs.begin(), pairs.end());

    std::vector<T> left_members;
    std::vector<T> right_members;

    left_members.reserve(n_pairs);
    right_members.reserve(n_pairs);

    for (const auto [a,b] : pairs) {
        left_members.emplace_back(a);
        right_members.emplace_back(b);
    }

    return std::make_pair(left_members, right_members);
}

// Maybe I could use c++20's std::popcount instead??
static inline int32_t popcnt128(__m128i n) {
    const __m128i n_hi = _mm_unpackhi_epi64(n, n);
    return __builtin_popcountll(_mm_cvtsi128_si64(n)) + __builtin_popcountll(_mm_cvtsi128_si64(n_hi));
}

static void BM_fp(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    std::vector<int8_t> phenotypes_ = make_random_bools();

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

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cncn);
}

BENCHMARK(BM_fp);

static void BM_nochecks_fp(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    std::vector<int8_t> phenotypes_ = make_random_bools();

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

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cncn);}

BENCHMARK(BM_nochecks_fp);

static void BM_int(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    std::vector<int8_t> phenotypes_ = make_random_bools();

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

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cncn);}

BENCHMARK(BM_int);

static void BM_nochecks_int(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    std::vector<int8_t> phenotypes_ = make_random_bools();

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

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cncn);}

BENCHMARK(BM_nochecks_int);

static void BM_nochecks_bool(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    auto pheno_bytes = make_random_bools();

    std::vector<bool> phenotypes_;
    phenotypes_.reserve(NUM_PHENOS);

    for (auto byte : pheno_bytes) {
        phenotypes_.emplace_back(byte);
    }

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    for (auto _ : state) {
        for (const auto p : pairs) {
            const auto [p1, p2] = p;
            const auto x = phenotypes_[p1];
            const auto y = phenotypes_[p2];

            cscs += x & y;
            cscn += x ^ y;
            cncn += !(x || y);
        }
    }

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cncn);}

BENCHMARK(BM_nochecks_bool);

static void BM_nochecks_valarray_bool(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    auto pheno_bytes = make_random_bools();

    std::valarray<bool> phenotypes_(NUM_PHENOS);

    for (size_t i = 0; i < pheno_bytes.size(); i++) {
        phenotypes_[i] = pheno_bytes[i];
    }

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    for (auto _ : state) {
        for (const auto p : pairs) {
            const auto [p1, p2] = p;
            const auto x = phenotypes_[p1];
            const auto y = phenotypes_[p2];

            cscs += x & y;
            cscn += x ^ y;
            cncn += !(x || y);
        }
    }

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cncn);}

BENCHMARK(BM_nochecks_valarray_bool);

static void BM_nochecks_bool_nocncn(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    std::vector<int8_t> phenotypes_ = make_random_bools();

    int64_t cscs = 0;
    int64_t cscn = 0;

    for (auto _ : state) {
        for (const auto p : pairs) {
            const auto [p1, p2] = p;
            const auto x = phenotypes_[p1];
            const auto y = phenotypes_[p2];

            cscs += x & y;
            cscn += x ^ y;
        }
    }

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
}

BENCHMARK(BM_nochecks_bool_nocncn);

static void BM_nochecks_bool_32_bit_pairs_nocncn(benchmark::State &state) {
    auto pairs = make_pairs<int32_t>();
    std::vector<int8_t> phenotypes_ = make_random_bools();

    int64_t cscs = 0;
    int64_t cscn = 0;

    for (auto _ : state) {
        for (const auto p : pairs) {
            const auto [p1, p2] = p;
            const auto x = phenotypes_[p1];
            const auto y = phenotypes_[p2];

            cscs += x & y;
            cscn += x ^ y;
        }
    }

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
}

BENCHMARK(BM_nochecks_bool_32_bit_pairs_nocncn);

static void BM_fastest_two_accumulators(benchmark::State &state) {
    auto pairs = make_pairs<int32_t>();
    std::vector<int8_t> phenotypes_ = make_random_bools();

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cscs2 = 0;
    int64_t cscn2 = 0;

    for (auto _ : state) {
        for (size_t i = 0; i < pairs.size() - 1; i += 2) {

            const auto [p1_a, p1_b] = pairs[i];
            const auto x1 = phenotypes_[p1_a];
            const auto y1 = phenotypes_[p1_b];

            cscs += x1 & y1;
            cscn += x1 ^ y1;

            const auto [p2_a, p2_b] = pairs[i + 1];
            auto x2 = phenotypes_[p2_a];
            auto y2 = phenotypes_[p2_b];

            cscs2 += x2 & y2;
            cscn2 += x2 ^ y2;
        }
    }

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
    benchmark::DoNotOptimize(cscs2);
    benchmark::DoNotOptimize(cscn2);
}

BENCHMARK(BM_fastest_two_accumulators);

static void BM_access_only(benchmark::State &state) {
    std::vector<std::pair<size_t, size_t>> pairs = make_pairs();
    std::vector<int8_t> phenotypes_ = make_random_bools();

    for (auto _ : state) {
        for (const auto p : pairs) {
            const auto [p1, p2] = p;
            const auto x = phenotypes_[p1];
            const auto y = phenotypes_[p2];

            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(y);
        }
    }
}

BENCHMARK(BM_access_only);

static void BM_vectorized_access_only(benchmark::State &state) {
    auto [left_members, right_members] = make_pairs_struct_of_arrays<int32_t>();
    std::vector<int8_t> phenotypes_ = make_random_bools();

    int64_t cscs = 0;
    int64_t cscn = 0;

    for (auto _ : state) {
        
        for (size_t i = 0; i < left_members.size() - 15; i+= 16) {
            auto left_addresses = _mm512_loadu_si512(&left_members[i]);
            auto right_addresses = _mm512_loadu_si512(&right_members[i]);

            auto lefts = _mm512_i32gather_epi32(left_addresses, phenotypes_.data(), 1);
            auto rights = _mm512_i32gather_epi32(right_addresses, phenotypes_.data(), 1);
            
            auto left_packed = _mm512_cvtepi32_epi8(lefts);
            auto right_packed = _mm512_cvtepi32_epi8(rights);

            benchmark::DoNotOptimize(left_packed);
            benchmark::DoNotOptimize(right_packed);
        }
    }
}

BENCHMARK(BM_vectorized_access_only);

static void BM_vectorized(benchmark::State &state) {
    auto [left_members, right_members] = make_pairs_struct_of_arrays<int32_t>();
    std::vector<int8_t> phenotypes_ = make_random_bools();

    int64_t cscs = 0;
    int64_t cscn = 0;

    for (auto _ : state) {
        for (size_t i = 0; i < left_members.size() - 15; i+= 16) {
            auto left_addresses = _mm512_loadu_si512(&left_members[i]);
            auto right_addresses = _mm512_loadu_si512(&right_members[i]);

            auto lefts = _mm512_i32gather_epi32(left_addresses, phenotypes_.data(), 1);
            auto rights = _mm512_i32gather_epi32(right_addresses, phenotypes_.data(), 1);
            
            auto left_packed = _mm512_cvtepi32_epi8(lefts);
            auto right_packed = _mm512_cvtepi32_epi8(rights);

            auto cscs_batch = _mm_and_si128(left_packed, right_packed);
            auto cscn_batch = _mm_xor_si128(left_packed, right_packed);

            cscs += popcnt128(cscs_batch);
            cscn += popcnt128(cscn_batch);
        }
    }

    benchmark::DoNotOptimize(cscs);
    benchmark::DoNotOptimize(cscn);
}

BENCHMARK(BM_vectorized);

// This could still be made faster by doing aligned loads rather than unaligned

BENCHMARK_MAIN();
