#ifndef CARVAIBD_TRANSPOSED_PHENOTYPES_HPP
#define CARVAIBD_TRANSPOSED_PHENOTYPES_HPP

#include "types.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <mutex>
#include <stdexcept>
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

    uint64_t* row(size_t sample) {
        return data.data() + sample * stride_words;
    }

    BitPackedPhenotypes(size_t n_samples_, size_t n_perms_)
        : n_samples(n_samples_), n_perms(n_perms_), n_total(n_perms_ + 1) {
        // Align stride to 8 words (512 bits) for AVX-512 friendly access.
        stride_words = align_up((n_total + 63) / 64, 8);
        data.resize(n_samples * stride_words, 0);
    }

    BitPackedPhenotypes(const std::vector<pheno_vector>& perms, size_t n_samples_)
        : BitPackedPhenotypes(n_samples_, perms.size() - 1) {

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

/**
 * Thread-safe destination for word-sized blocks of bit-packed permutations.
 * A lock is taken once per completed word, not once per phenotype value.
 */
class BitPackedPermutationWriter {
    BitPackedPhenotypes& output;
    std::vector<std::mutex> word_mutexes;

public:
    explicit BitPackedPermutationWriter(BitPackedPhenotypes& output_)
        : output(output_), word_mutexes((output_.n_total + 63) / 64) {}

    void commit(size_t word, const std::vector<uint64_t>& values) {
        if (word >= word_mutexes.size() || values.size() != output.n_samples) {
            throw std::out_of_range("Invalid bit-packed permutation word");
        }
        std::lock_guard<std::mutex> lock(word_mutexes[word]);
        for (size_t sample = 0; sample < output.n_samples; ++sample) {
            output.row(sample)[word] |= values[sample];
        }
    }
};

/**
 * Accumulates a contiguous range of permutations into 64-bit sample words.
 * Multiple buffers may safely commit disjoint bits of the same word.
 */
class BitPackedPermutationBuffer {
    BitPackedPermutationWriter& writer;
    size_t n_samples;
    size_t n_total;
    size_t current_word = 0;
    size_t current_permutation = 0;
    bool active = false;
    std::vector<uint64_t> values;

    void flush() {
        if (active) {
            writer.commit(current_word, values);
            std::fill(values.begin(), values.end(), uint64_t{0});
        }
    }

public:
    BitPackedPermutationBuffer(BitPackedPermutationWriter& writer_,
                               const BitPackedPhenotypes& output)
        : writer(writer_), n_samples(output.n_samples), n_total(output.n_total),
          values(n_samples, 0) {}

    void begin_permutation(size_t permutation) {
        if (permutation >= n_total) {
            throw std::out_of_range("Permutation index exceeds bit-packed storage");
        }
        const size_t word = permutation / 64;
        if (active && word != current_word) {
            flush();
            active = false;
        }
        current_word = word;
        current_permutation = permutation;
        active = true;
    }

    void set_case(size_t sample) {
        if (!active || sample >= n_samples) {
            throw std::out_of_range("Invalid bit-packed phenotype sample");
        }
        values[sample] |= uint64_t{1} << (current_permutation % 64);
    }

    template<typename T>
    void append(size_t permutation, const T& phenotype) {
        if (phenotype.size() != n_samples) {
            throw std::runtime_error("Permutation has the wrong sample count");
        }
        begin_permutation(permutation);
        for (size_t sample = 0; sample < n_samples; ++sample) {
            if (phenotype[sample]) {
                set_case(sample);
            }
        }
    }

    void finish() {
        flush();
        active = false;
    }
};

#endif // CARVAIBD_TRANSPOSED_PHENOTYPES_HPP
