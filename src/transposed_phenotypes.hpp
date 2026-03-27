#ifndef CARVAIBD_TRANSPOSED_PHENOTYPES_HPP
#define CARVAIBD_TRANSPOSED_PHENOTYPES_HPP

#include "types.hpp"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

inline size_t align_up(size_t n, size_t alignment) {
    return (n + alignment - 1) & ~(alignment - 1);
}

struct TransposedPhenotypes {
    std::vector<uint8_t> data;
    size_t stride;     // bytes per sample row, aligned to 64
    size_t n_samples;
    size_t n_perms;    // number of permutations (not counting original)

    const uint8_t* row(size_t sample) const {
        return data.data() + sample * stride;
    }

    TransposedPhenotypes(const std::vector<pheno_vector>& perms, size_t n_samples_)
        : n_samples(n_samples_), n_perms(perms.size() - 1) {
        size_t n_cols = perms.size(); // nperms + 1 (original at index 0)
        stride = align_up(n_cols, 64);
        data.resize(n_samples * stride, 0);

        for (size_t perm = 0; perm < n_cols; ++perm) {
            const auto& pv = perms[perm];
            for (size_t s = 0; s < n_samples; ++s) {
                data[s * stride + perm] = static_cast<uint8_t>(pv[s]);
            }
        }
    }
};

/**
 * Bit-packed transposed phenotype storage.
 * Each sample's permutation values are packed as bits into uint64_t words.
 * Layout: data[sample * stride_words + word] contains 64 permutation bits.
 * Bit k of word w corresponds to permutation (w * 64 + k).
 */
struct BitPackedPhenotypes {
    std::vector<uint64_t> data;
    size_t stride_words;   // uint64_t words per sample row, aligned to 8 (512 bits)
    size_t n_samples;
    size_t n_perms;        // number of permutations (not counting original)
    size_t n_total;        // nperms + 1

    const uint64_t* row(size_t sample) const {
        return data.data() + sample * stride_words;
    }

    BitPackedPhenotypes(const std::vector<pheno_vector>& perms, size_t n_samples_)
        : n_samples(n_samples_), n_perms(perms.size() - 1), n_total(perms.size()) {
        // Align stride to 8 words (512 bits) for AVX-512 friendly access
        stride_words = align_up((n_total + 63) / 64, 8);
        data.resize(n_samples * stride_words, 0);

        for (size_t perm = 0; perm < n_total; ++perm) {
            size_t word = perm / 64;
            size_t bit  = perm % 64;
            uint64_t mask = uint64_t(1) << bit;
            const auto& pv = perms[perm];
            for (size_t s = 0; s < n_samples; ++s) {
                if (pv[s]) {
                    data[s * stride_words + word] |= mask;
                }
            }
        }
    }
};

#endif // CARVAIBD_TRANSPOSED_PHENOTYPES_HPP
