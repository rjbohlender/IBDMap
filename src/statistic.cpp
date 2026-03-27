//
// Created by Bohlender,Ryan James on 9/23/19.
//

#include "statistic.hpp"
#include "split.hpp"
#include <cstring>
#include <fmt/include/fmt/ostream.h>
#include <utility>

#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif

template<typename T>
Statistic<T>::Statistic(arma::SpCol<int32_t> data_,
                        Breakpoint bp_,
                        std::shared_ptr<Indexer<T>> indexer_,
                        std::shared_ptr<Reporter> reporter_,
                        uint64_t seq_,
                        Parameters params_,
                        std::shared_ptr<std::vector<T>> phenotypes_,
                        std::shared_ptr<TransposedPhenotypes> transposed_,
                        std::shared_ptr<BitPackedPhenotypes> bitpacked_) : data(std::move(data_)), indexer(std::move(indexer_)),
                                                                       seq(seq_),
                                                                       params(std::move(params_)),
                                                                       bp(std::move(bp_)), reporter(std::move(reporter_)),
                                                                       phenotypes(std::move(phenotypes_)),
                                                                       transposed(std::move(transposed_)),
                                                                       bitpacked(std::move(bitpacked_)) {
    if (params.enable_testing) {
        test_statistic();
    }
}

template<typename T>
double
Statistic<T>::calculate(T &phenotypes_, bool original_) {
    double statistic = 0;

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    if (pairs.first.empty()) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            auto [left, right] = (*indexer).back_translate(it.row());
            pairs.first.emplace_back(left);
            pairs.second.emplace_back(right);
        }
    }

    const size_t n = pairs.first.size();

#if defined(__AVX512F__) && defined(__AVX512BW__)
    if constexpr (std::is_same<T, pheno_vector>::value) {
        // Pre-materialize phenotype values into contiguous arrays
        if (left_phenos.size() != n) {
            left_phenos.resize(n + 63);  // pad for SIMD overread
            right_phenos.resize(n + 63);
        }
        for (size_t i = 0; i < n; ++i) {
            left_phenos[i] = static_cast<uint8_t>(phenotypes_[pairs.first[i]]);
            right_phenos[i] = static_cast<uint8_t>(phenotypes_[pairs.second[i]]);
        }

        static bool printed = false;
        if (!printed) {
            fmt::print(std::cerr, "Using AVX-512 path for statistic calculation\n");
            printed = true;
        }

        __m512i cscs_acc = _mm512_setzero_si512();
        __m512i cscn_acc = _mm512_setzero_si512();
        size_t i = 0;
        size_t batch = 0;

        for (; i + 63 < n; i += 64) {
            __m512i l = _mm512_loadu_si512(left_phenos.data() + i);
            __m512i r = _mm512_loadu_si512(right_phenos.data() + i);
            cscs_acc = _mm512_add_epi8(cscs_acc, _mm512_and_si512(l, r));
            cscn_acc = _mm512_add_epi8(cscn_acc, _mm512_xor_si512(l, r));
            batch++;
            if (batch == 255) {
                // Flush byte accumulators using sad_epu8 (sum of absolute differences with zero)
                __m512i zero = _mm512_setzero_si512();
                __m512i cscs_sad = _mm512_sad_epu8(cscs_acc, zero);  // 8 × uint64 partial sums
                __m512i cscn_sad = _mm512_sad_epu8(cscn_acc, zero);
                // Reduce 8 uint64 lanes
                cscs += _mm512_reduce_add_epi64(cscs_sad);
                cscn += _mm512_reduce_add_epi64(cscn_sad);
                cscs_acc = _mm512_setzero_si512();
                cscn_acc = _mm512_setzero_si512();
                batch = 0;
            }
        }
        // Final flush
        if (batch > 0) {
            __m512i zero = _mm512_setzero_si512();
            cscs += _mm512_reduce_add_epi64(_mm512_sad_epu8(cscs_acc, zero));
            cscn += _mm512_reduce_add_epi64(_mm512_sad_epu8(cscn_acc, zero));
        }
        // Scalar tail
        for (; i < n; ++i) {
            cscs += left_phenos[i] & right_phenos[i];
            cscn += left_phenos[i] ^ right_phenos[i];
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            const auto x = phenotypes_[pairs.first[i]];
            const auto y = phenotypes_[pairs.second[i]];
            cscs += x & y;
            cscn += x ^ y;
        }
    }

#elif defined(__AVX2__)
    if constexpr (std::is_same<T, pheno_vector>::value) {
        // Pre-materialize phenotype values into contiguous arrays
        if (left_phenos.size() != n) {
            left_phenos.resize(n + 31);  // pad for SIMD overread
            right_phenos.resize(n + 31);
        }
        for (size_t i = 0; i < n; ++i) {
            left_phenos[i] = static_cast<uint8_t>(phenotypes_[pairs.first[i]]);
            right_phenos[i] = static_cast<uint8_t>(phenotypes_[pairs.second[i]]);
        }

        static bool printed = false;
        if (!printed) {
            fmt::print(std::cerr, "Using AVX2 path for statistic calculation\n");
            printed = true;
        }

        __m256i cscs_acc = _mm256_setzero_si256();
        __m256i cscn_acc = _mm256_setzero_si256();
        __m256i zero = _mm256_setzero_si256();
        size_t i = 0;
        size_t batch = 0;

        for (; i + 31 < n; i += 32) {
            __m256i l = _mm256_loadu_si256((const __m256i*)(left_phenos.data() + i));
            __m256i r = _mm256_loadu_si256((const __m256i*)(right_phenos.data() + i));
            cscs_acc = _mm256_add_epi8(cscs_acc, _mm256_and_si256(l, r));
            cscn_acc = _mm256_add_epi8(cscn_acc, _mm256_xor_si256(l, r));
            batch++;
            if (batch == 255) {
                // Flush: sad_epu8 sums 8 bytes per 64-bit lane → 4 partial sums
                __m256i cscs_sad = _mm256_sad_epu8(cscs_acc, zero);
                __m256i cscn_sad = _mm256_sad_epu8(cscn_acc, zero);
                // Extract and accumulate 4 uint64 lanes
                cscs += _mm256_extract_epi64(cscs_sad, 0) + _mm256_extract_epi64(cscs_sad, 1)
                      + _mm256_extract_epi64(cscs_sad, 2) + _mm256_extract_epi64(cscs_sad, 3);
                cscn += _mm256_extract_epi64(cscn_sad, 0) + _mm256_extract_epi64(cscn_sad, 1)
                      + _mm256_extract_epi64(cscn_sad, 2) + _mm256_extract_epi64(cscn_sad, 3);
                cscs_acc = _mm256_setzero_si256();
                cscn_acc = _mm256_setzero_si256();
                batch = 0;
            }
        }
        // Final flush
        if (batch > 0) {
            __m256i cscs_sad = _mm256_sad_epu8(cscs_acc, zero);
            __m256i cscn_sad = _mm256_sad_epu8(cscn_acc, zero);
            cscs += _mm256_extract_epi64(cscs_sad, 0) + _mm256_extract_epi64(cscs_sad, 1)
                  + _mm256_extract_epi64(cscs_sad, 2) + _mm256_extract_epi64(cscs_sad, 3);
            cscn += _mm256_extract_epi64(cscn_sad, 0) + _mm256_extract_epi64(cscn_sad, 1)
                  + _mm256_extract_epi64(cscn_sad, 2) + _mm256_extract_epi64(cscn_sad, 3);
        }
        // Scalar tail
        for (; i < n; ++i) {
            cscs += left_phenos[i] & right_phenos[i];
            cscn += left_phenos[i] ^ right_phenos[i];
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            const auto x = phenotypes_[pairs.first[i]];
            const auto y = phenotypes_[pairs.second[i]];
            cscs += x & y;
            cscn += x ^ y;
        }
    }

#else
    for (size_t i = 0; i < n; ++i) {
        const auto x = phenotypes_[pairs.first[i]];
        const auto y = phenotypes_[pairs.second[i]];
        cscs += x & y;
        cscn += x ^ y;
    }
#endif

    cncn = pairs.first.size() - cscs - cscn;

    if (original_) {
        orig_cscs = static_cast<double>(cscs);
        orig_cscn = static_cast<double>(cscn);
        orig_cncn = static_cast<double>(cncn);
    }

    if (params.print_debug) {
        reporter->submit(seq, fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                     bp.breakpoint.first,
                                     bp.breakpoint.second,
                                     (*indexer).case_case,
                                     (*indexer).case_cont,
                                     (*indexer).cont_cont,
                                     cscs,
                                     cscn,
                                     cncn), true);
    }

    // I have some numerical concerns here. This may be prone to catastrophic cancellation and could result in incorrect results.
    // Going to transform it so that we can do everything with integers until we have a single division
    if (params.contcont) {
        statistic = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale - cncn * lcm_cncn_scale) / lcm_common;
    } else if (params.cscs_only) {
        statistic = static_cast<double>(cscs) / static_cast<double>((*indexer).case_case);
    } else {
        statistic = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale) / lcm_common;
    }

    return statistic;
}

template<typename T>
void Statistic<T>::run() {
    setup_lcm();
    if (bitpacked) {
        permute_bulk_bitpacked();
    } else if (transposed) {
        permute_bulk();
    } else {
        initialize();
        permute();
    }

    // Output
    std::stringstream ss;
    build_output(ss);
    if (params.verbose) {
        fmt::print(std::cerr, "Finished {}\t{}\n", bp.breakpoint.first, bp.breakpoint.second);
    }

    reporter->submit(seq, ss.str(), false);
    cleanup();
    done = true;
}

template<typename T>
void Statistic<T>::build_output(std::stringstream &ss) {
    double cscs = (*indexer).case_case;
    double cscn = (*indexer).case_cont;
    double cncn = (*indexer).cont_cont;
    fmt::print(ss, "{}\t{}\t", bp.breakpoint.first, bp.breakpoint.second);
    fmt::print(ss, "{}\t{}\t{}\t{}", orig_cscs / cscs, orig_cscn / cscn, orig_cncn / cncn, original);
    for (int i = 0; i < params.nperms; i++) {
        fmt::print(ss, "\t{}", permuted[i]);
    }
    fmt::print(ss, "\n");
}

template<typename T>
void Statistic<T>::permute() {
    permuted.reserve(params.nperms);
    double val;
    for (int i = 1; i <= params.nperms; i++) {
        val = calculate(phenotypes->at(i), false);
        permuted.push_back(val);
    }
}

template<typename T>
void Statistic<T>::permute_bulk() {
    if (pairs.first.empty()) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            auto [left, right] = (*indexer).back_translate(it.row());
            pairs.first.emplace_back(left);
            pairs.second.emplace_back(right);
        }
    }

    const size_t n_perms = transposed->n_perms;
    const size_t stride = transposed->stride;
    const size_t n_pairs = pairs.first.size();

    // Per-permutation accumulators: index 0 = original, 1..n_perms = permutations
    // We use int32_t since pair counts fit in 32 bits
    std::vector<int32_t> cscs_counts(stride, 0);
    std::vector<int32_t> cscn_counts(stride, 0);

    // Process pairs in batches of 255 to use byte accumulators without overflow
    constexpr size_t BATCH = 255;

#if defined(__AVX512F__) && defined(__AVX512BW__)
    static bool simd_printed = false;
    if (!simd_printed) {
        fmt::print(std::cerr, "Using AVX-512 bulk permutation path\n");
        simd_printed = true;
    }
    constexpr size_t VEC_WIDTH = 64;
    size_t n_vecs = stride / VEC_WIDTH;
    std::vector<__m512i> cscs_bytes(n_vecs);
    std::vector<__m512i> cscn_bytes(n_vecs);

    for (size_t batch_start = 0; batch_start < n_pairs; batch_start += BATCH) {
        size_t batch_end = std::min(batch_start + BATCH, n_pairs);

        for (size_t v = 0; v < n_vecs; ++v) {
            cscs_bytes[v] = _mm512_setzero_si512();
            cscn_bytes[v] = _mm512_setzero_si512();
        }

        for (size_t pi = batch_start; pi < batch_end; ++pi) {
            const uint8_t* L = transposed->row(pairs.first[pi]);
            const uint8_t* R = transposed->row(pairs.second[pi]);
            for (size_t v = 0; v < n_vecs; ++v) {
                __m512i l = _mm512_loadu_si512(L + v * VEC_WIDTH);
                __m512i r = _mm512_loadu_si512(R + v * VEC_WIDTH);
                cscs_bytes[v] = _mm512_add_epi8(cscs_bytes[v], _mm512_and_si512(l, r));
                cscn_bytes[v] = _mm512_add_epi8(cscn_bytes[v], _mm512_xor_si512(l, r));
            }
        }

        // Flush byte accumulators to 32-bit counts
        for (size_t v = 0; v < n_vecs; ++v) {
            alignas(64) uint8_t cscs_tmp[64], cscn_tmp[64];
            _mm512_store_si512(cscs_tmp, cscs_bytes[v]);
            _mm512_store_si512(cscn_tmp, cscn_bytes[v]);
            size_t base = v * VEC_WIDTH;
            for (size_t b = 0; b < VEC_WIDTH; ++b) {
                cscs_counts[base + b] += cscs_tmp[b];
                cscn_counts[base + b] += cscn_tmp[b];
            }
        }
    }

#elif defined(__AVX2__)
    static bool simd_printed = false;
    if (!simd_printed) {
        fmt::print(std::cerr, "Using AVX2 bulk permutation path\n");
        simd_printed = true;
    }
    constexpr size_t VEC_WIDTH = 32;
    size_t n_vecs = stride / VEC_WIDTH;
    std::vector<__m256i> cscs_bytes(n_vecs);
    std::vector<__m256i> cscn_bytes(n_vecs);

    for (size_t batch_start = 0; batch_start < n_pairs; batch_start += BATCH) {
        size_t batch_end = std::min(batch_start + BATCH, n_pairs);

        for (size_t v = 0; v < n_vecs; ++v) {
            cscs_bytes[v] = _mm256_setzero_si256();
            cscn_bytes[v] = _mm256_setzero_si256();
        }

        for (size_t pi = batch_start; pi < batch_end; ++pi) {
            const uint8_t* L = transposed->row(pairs.first[pi]);
            const uint8_t* R = transposed->row(pairs.second[pi]);
            for (size_t v = 0; v < n_vecs; ++v) {
                __m256i l = _mm256_loadu_si256((const __m256i*)(L + v * VEC_WIDTH));
                __m256i r = _mm256_loadu_si256((const __m256i*)(R + v * VEC_WIDTH));
                cscs_bytes[v] = _mm256_add_epi8(cscs_bytes[v], _mm256_and_si256(l, r));
                cscn_bytes[v] = _mm256_add_epi8(cscn_bytes[v], _mm256_xor_si256(l, r));
            }
        }

        // Flush byte accumulators to 32-bit counts
        for (size_t v = 0; v < n_vecs; ++v) {
            alignas(32) uint8_t cscs_tmp[32], cscn_tmp[32];
            _mm256_store_si256((__m256i*)cscs_tmp, cscs_bytes[v]);
            _mm256_store_si256((__m256i*)cscn_tmp, cscn_bytes[v]);
            size_t base = v * VEC_WIDTH;
            for (size_t b = 0; b < VEC_WIDTH; ++b) {
                cscs_counts[base + b] += cscs_tmp[b];
                cscn_counts[base + b] += cscn_tmp[b];
            }
        }
    }

#else
    // Scalar fallback — no byte overflow concern with int32_t accumulators
    for (size_t pi = 0; pi < n_pairs; ++pi) {
        const uint8_t* L = transposed->row(pairs.first[pi]);
        const uint8_t* R = transposed->row(pairs.second[pi]);
        for (size_t j = 0; j < n_perms + 1; ++j) {
            cscs_counts[j] += L[j] & R[j];
            cscn_counts[j] += L[j] ^ R[j];
        }
    }
#endif

    // Extract original statistic (index 0)
    {
        int64_t cscs = cscs_counts[0];
        int64_t cscn = cscn_counts[0];
        int64_t cncn = static_cast<int64_t>(n_pairs) - cscs - cscn;
        orig_cscs = static_cast<double>(cscs);
        orig_cscn = static_cast<double>(cscn);
        orig_cncn = static_cast<double>(cncn);

        if (params.print_debug) {
            reporter->submit(seq, fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                         bp.breakpoint.first,
                                         bp.breakpoint.second,
                                         (*indexer).case_case,
                                         (*indexer).case_cont,
                                         (*indexer).cont_cont,
                                         cscs, cscn, cncn), true);
        }

        if (params.contcont) {
            original = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale - cncn * lcm_cncn_scale) / lcm_common;
        } else if (params.cscs_only) {
            original = static_cast<double>(cscs) / static_cast<double>((*indexer).case_case);
        } else {
            original = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale) / lcm_common;
        }
    }

    // Convert per-permutation counts to statistics
    permuted.resize(n_perms);
    for (size_t i = 0; i < n_perms; ++i) {
        int64_t cscs = cscs_counts[i + 1];
        int64_t cscn = cscn_counts[i + 1];
        if (params.contcont) {
            int64_t cncn = static_cast<int64_t>(n_pairs) - cscs - cscn;
            permuted[i] = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale - cncn * lcm_cncn_scale) / lcm_common;
        } else if (params.cscs_only) {
            permuted[i] = static_cast<double>(cscs) / static_cast<double>((*indexer).case_case);
        } else {
            permuted[i] = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale) / lcm_common;
        }
    }
}

template<typename T>
void Statistic<T>::permute_bulk_bitpacked() {
    if (pairs.first.empty()) {
        for (auto it = data.begin(); it != data.end(); ++it) {
            auto [left, right] = (*indexer).back_translate(it.row());
            pairs.first.emplace_back(left);
            pairs.second.emplace_back(right);
        }
    }

    const size_t n_total = bitpacked->n_total;
    const size_t n_perms = bitpacked->n_perms;
    const size_t stride_words = bitpacked->stride_words;
    const size_t n_pairs = pairs.first.size();

    std::vector<int32_t> cscs_counts(n_total, 0);
    std::vector<int32_t> cscn_counts(n_total, 0);

    constexpr size_t BATCH = 255;
    constexpr int N_PLANES = 8;

#if defined(__AVX512F__)
    static bool simd_printed = false;
    if (!simd_printed) {
        fmt::print(std::cerr, "Using AVX-512 bit-packed permutation path\n");
        simd_printed = true;
    }

    constexpr size_t WORDS_PER_VEC = 8;
    const size_t n_vecs = stride_words / WORDS_PER_VEC;

    for (size_t vw = 0; vw < n_vecs; ++vw) {
        __m512i cscs_planes[N_PLANES];
        __m512i cscn_planes[N_PLANES];
        for (int k = 0; k < N_PLANES; ++k) {
            cscs_planes[k] = _mm512_setzero_si512();
            cscn_planes[k] = _mm512_setzero_si512();
        }

        size_t batch_count = 0;
        for (size_t pi = 0; pi < n_pairs; ++pi) {
            const uint64_t* L = bitpacked->row(pairs.first[pi]) + vw * WORDS_PER_VEC;
            const uint64_t* R = bitpacked->row(pairs.second[pi]) + vw * WORDS_PER_VEC;

            __m512i lv = _mm512_loadu_si512(L);
            __m512i rv = _mm512_loadu_si512(R);
            __m512i both = _mm512_and_si512(lv, rv);
            __m512i diff = _mm512_xor_si512(lv, rv);

            for (int k = 0; k < N_PLANES; ++k) {
                __m512i new_plane = _mm512_xor_si512(cscs_planes[k], both);
                both = _mm512_and_si512(cscs_planes[k], both);
                cscs_planes[k] = new_plane;
            }
            for (int k = 0; k < N_PLANES; ++k) {
                __m512i new_plane = _mm512_xor_si512(cscn_planes[k], diff);
                diff = _mm512_and_si512(cscn_planes[k], diff);
                cscn_planes[k] = new_plane;
            }

            batch_count++;
            if (batch_count == BATCH) {
                // Flush planes to int32 counts
                for (int k = 0; k < N_PLANES; ++k) {
                    alignas(64) uint64_t cscs_buf[WORDS_PER_VEC], cscn_buf[WORDS_PER_VEC];
                    _mm512_store_si512(cscs_buf, cscs_planes[k]);
                    _mm512_store_si512(cscn_buf, cscn_planes[k]);
                    for (size_t elem = 0; elem < WORDS_PER_VEC; ++elem) {
                        size_t base = (vw * WORDS_PER_VEC + elem) * 64;
                        size_t max_bit = (base < n_total) ? std::min(size_t(64), n_total - base) : 0;
                        uint64_t cscs_word = cscs_buf[elem];
                        uint64_t cscn_word = cscn_buf[elem];
                        for (size_t b = 0; b < max_bit; ++b) {
                            cscs_counts[base + b] += static_cast<int32_t>((cscs_word >> b) & 1) << k;
                            cscn_counts[base + b] += static_cast<int32_t>((cscn_word >> b) & 1) << k;
                        }
                    }
                    cscs_planes[k] = _mm512_setzero_si512();
                    cscn_planes[k] = _mm512_setzero_si512();
                }
                batch_count = 0;
            }
        }

        // Final flush
        if (batch_count > 0) {
            for (int k = 0; k < N_PLANES; ++k) {
                alignas(64) uint64_t cscs_buf[WORDS_PER_VEC], cscn_buf[WORDS_PER_VEC];
                _mm512_store_si512(cscs_buf, cscs_planes[k]);
                _mm512_store_si512(cscn_buf, cscn_planes[k]);
                for (size_t elem = 0; elem < WORDS_PER_VEC; ++elem) {
                    size_t base = (vw * WORDS_PER_VEC + elem) * 64;
                    size_t max_bit = (base < n_total) ? std::min(size_t(64), n_total - base) : 0;
                    uint64_t cscs_word = cscs_buf[elem];
                    uint64_t cscn_word = cscn_buf[elem];
                    for (size_t b = 0; b < max_bit; ++b) {
                        cscs_counts[base + b] += static_cast<int32_t>((cscs_word >> b) & 1) << k;
                        cscn_counts[base + b] += static_cast<int32_t>((cscn_word >> b) & 1) << k;
                    }
                }
            }
        }
    }

#elif defined(__AVX2__)
    static bool simd_printed = false;
    if (!simd_printed) {
        fmt::print(std::cerr, "Using AVX2 bit-packed permutation path\n");
        simd_printed = true;
    }

    constexpr size_t WORDS_PER_VEC = 4;
    const size_t n_vecs = stride_words / WORDS_PER_VEC;

    for (size_t vw = 0; vw < n_vecs; ++vw) {
        __m256i cscs_planes[N_PLANES];
        __m256i cscn_planes[N_PLANES];
        for (int k = 0; k < N_PLANES; ++k) {
            cscs_planes[k] = _mm256_setzero_si256();
            cscn_planes[k] = _mm256_setzero_si256();
        }

        size_t batch_count = 0;
        for (size_t pi = 0; pi < n_pairs; ++pi) {
            const uint64_t* L = bitpacked->row(pairs.first[pi]) + vw * WORDS_PER_VEC;
            const uint64_t* R = bitpacked->row(pairs.second[pi]) + vw * WORDS_PER_VEC;

            __m256i lv = _mm256_loadu_si256((const __m256i*)L);
            __m256i rv = _mm256_loadu_si256((const __m256i*)R);
            __m256i both = _mm256_and_si256(lv, rv);
            __m256i diff = _mm256_xor_si256(lv, rv);

            for (int k = 0; k < N_PLANES; ++k) {
                __m256i new_plane = _mm256_xor_si256(cscs_planes[k], both);
                both = _mm256_and_si256(cscs_planes[k], both);
                cscs_planes[k] = new_plane;
            }
            for (int k = 0; k < N_PLANES; ++k) {
                __m256i new_plane = _mm256_xor_si256(cscn_planes[k], diff);
                diff = _mm256_and_si256(cscn_planes[k], diff);
                cscn_planes[k] = new_plane;
            }

            batch_count++;
            if (batch_count == BATCH) {
                for (int k = 0; k < N_PLANES; ++k) {
                    alignas(32) uint64_t cscs_buf[WORDS_PER_VEC], cscn_buf[WORDS_PER_VEC];
                    _mm256_store_si256((__m256i*)cscs_buf, cscs_planes[k]);
                    _mm256_store_si256((__m256i*)cscn_buf, cscn_planes[k]);
                    for (size_t elem = 0; elem < WORDS_PER_VEC; ++elem) {
                        size_t base = (vw * WORDS_PER_VEC + elem) * 64;
                        size_t max_bit = (base < n_total) ? std::min(size_t(64), n_total - base) : 0;
                        uint64_t cscs_word = cscs_buf[elem];
                        uint64_t cscn_word = cscn_buf[elem];
                        for (size_t b = 0; b < max_bit; ++b) {
                            cscs_counts[base + b] += static_cast<int32_t>((cscs_word >> b) & 1) << k;
                            cscn_counts[base + b] += static_cast<int32_t>((cscn_word >> b) & 1) << k;
                        }
                    }
                    cscs_planes[k] = _mm256_setzero_si256();
                    cscn_planes[k] = _mm256_setzero_si256();
                }
                batch_count = 0;
            }
        }

        if (batch_count > 0) {
            for (int k = 0; k < N_PLANES; ++k) {
                alignas(32) uint64_t cscs_buf[WORDS_PER_VEC], cscn_buf[WORDS_PER_VEC];
                _mm256_store_si256((__m256i*)cscs_buf, cscs_planes[k]);
                _mm256_store_si256((__m256i*)cscn_buf, cscn_planes[k]);
                for (size_t elem = 0; elem < WORDS_PER_VEC; ++elem) {
                    size_t base = (vw * WORDS_PER_VEC + elem) * 64;
                    size_t max_bit = (base < n_total) ? std::min(size_t(64), n_total - base) : 0;
                    uint64_t cscs_word = cscs_buf[elem];
                    uint64_t cscn_word = cscn_buf[elem];
                    for (size_t b = 0; b < max_bit; ++b) {
                        cscs_counts[base + b] += static_cast<int32_t>((cscs_word >> b) & 1) << k;
                        cscn_counts[base + b] += static_cast<int32_t>((cscn_word >> b) & 1) << k;
                    }
                }
            }
        }
    }

#else
    // Scalar fallback: process one word (64 permutations) at a time
    for (size_t w = 0; w < stride_words; ++w) {
        uint64_t cscs_planes[N_PLANES] = {};
        uint64_t cscn_planes[N_PLANES] = {};

        size_t batch_count = 0;
        for (size_t pi = 0; pi < n_pairs; ++pi) {
            uint64_t both = bitpacked->row(pairs.first[pi])[w] & bitpacked->row(pairs.second[pi])[w];
            uint64_t diff = bitpacked->row(pairs.first[pi])[w] ^ bitpacked->row(pairs.second[pi])[w];

            for (int k = 0; k < N_PLANES; ++k) {
                uint64_t new_plane = cscs_planes[k] ^ both;
                both = cscs_planes[k] & both;
                cscs_planes[k] = new_plane;
                if (both == 0) break;
            }
            for (int k = 0; k < N_PLANES; ++k) {
                uint64_t new_plane = cscn_planes[k] ^ diff;
                diff = cscn_planes[k] & diff;
                cscn_planes[k] = new_plane;
                if (diff == 0) break;
            }

            batch_count++;
            if (batch_count == BATCH) {
                size_t base = w * 64;
                size_t max_bit = (base < n_total) ? std::min(size_t(64), n_total - base) : 0;
                for (size_t b = 0; b < max_bit; ++b) {
                    for (int k = 0; k < N_PLANES; ++k) {
                        cscs_counts[base + b] += static_cast<int32_t>((cscs_planes[k] >> b) & 1) << k;
                        cscn_counts[base + b] += static_cast<int32_t>((cscn_planes[k] >> b) & 1) << k;
                    }
                }
                std::memset(cscs_planes, 0, sizeof(cscs_planes));
                std::memset(cscn_planes, 0, sizeof(cscn_planes));
                batch_count = 0;
            }
        }

        if (batch_count > 0) {
            size_t base = w * 64;
            size_t max_bit = (base < n_total) ? std::min(size_t(64), n_total - base) : 0;
            for (size_t b = 0; b < max_bit; ++b) {
                for (int k = 0; k < N_PLANES; ++k) {
                    cscs_counts[base + b] += static_cast<int32_t>((cscs_planes[k] >> b) & 1) << k;
                    cscn_counts[base + b] += static_cast<int32_t>((cscn_planes[k] >> b) & 1) << k;
                }
            }
        }
    }
#endif

    // Extract original statistic (index 0)
    {
        int64_t cscs = cscs_counts[0];
        int64_t cscn = cscn_counts[0];
        int64_t cncn = static_cast<int64_t>(n_pairs) - cscs - cscn;
        orig_cscs = static_cast<double>(cscs);
        orig_cscn = static_cast<double>(cscn);
        orig_cncn = static_cast<double>(cncn);

        if (params.print_debug) {
            reporter->submit(seq, fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                         bp.breakpoint.first,
                                         bp.breakpoint.second,
                                         (*indexer).case_case,
                                         (*indexer).case_cont,
                                         (*indexer).cont_cont,
                                         cscs, cscn, cncn), true);
        }

        if (params.contcont) {
            original = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale - cncn * lcm_cncn_scale) / lcm_common;
        } else if (params.cscs_only) {
            original = static_cast<double>(cscs) / static_cast<double>((*indexer).case_case);
        } else {
            original = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale) / lcm_common;
        }
    }

    // Convert per-permutation counts to statistics
    permuted.resize(n_perms);
    for (size_t i = 0; i < n_perms; ++i) {
        int64_t cscs = cscs_counts[i + 1];
        int64_t cscn = cscn_counts[i + 1];
        if (params.contcont) {
            int64_t cncn = static_cast<int64_t>(n_pairs) - cscs - cscn;
            permuted[i] = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale - cncn * lcm_cncn_scale) / lcm_common;
        } else if (params.cscs_only) {
            permuted[i] = static_cast<double>(cscs) / static_cast<double>((*indexer).case_case);
        } else {
            permuted[i] = static_cast<double>(cscs * lcm_cscs_scale - cscn * lcm_cscn_scale) / lcm_common;
        }
    }
}

template<typename T>
void Statistic<T>::setup_lcm() {
    if (params.contcont) {
        lcm_common = static_cast<int64_t>(std::lcm((*indexer).case_case, std::lcm((*indexer).case_cont, (*indexer).cont_cont)));
        lcm_cscs_scale = lcm_common / static_cast<int64_t>((*indexer).case_case);
        lcm_cscn_scale = lcm_common / static_cast<int64_t>((*indexer).case_cont);
        lcm_cncn_scale = lcm_common / static_cast<int64_t>((*indexer).cont_cont);
    } else if (!params.cscs_only) {
        lcm_common = static_cast<int64_t>(std::lcm((*indexer).case_case, (*indexer).case_cont));
        lcm_cscs_scale = lcm_common / static_cast<int64_t>((*indexer).case_case);
        lcm_cscn_scale = lcm_common / static_cast<int64_t>((*indexer).case_cont);
    }
}

template<typename T>
void Statistic<T>::initialize() {
    original = calculate(phenotypes->at(0), true);
}

template<typename T>
void Statistic<T>::cleanup() {
    data.reset();
    permuted.clear();
    permuted.shrink_to_fit();
    left_phenos.clear();
    left_phenos.shrink_to_fit();
    right_phenos.clear();
    right_phenos.shrink_to_fit();
}

template<typename T>
void Statistic<T>::test_statistic() {
    // Test case setup
    std::vector<std::string> tsamples{
            "case1", "case2", "case3", "case4", "case5",
            "case6", "case7", "case8", "case9", "case10",
            "control1", "control2", "control3", "control4", "control5",
            "control6", "control7", "control8", "control9", "control10"};
    T tphenotypes_{
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Indexer tindexer(10, 10, tsamples, tphenotypes_);

    arma::sp_vec tdata(tsamples.size() * (tsamples.size() - 1.) / 2.);

    // 4 case-case pairs, 4 control-control pairs, 4 case-control pairs
    tdata(tindexer.translate("case1", "case2")) = 1.;
    tdata(tindexer.translate("case1", "case3")) = 1.;
    tdata(tindexer.translate("case1", "case4")) = 1.;
    tdata(tindexer.translate("case5", "case6")) = 1.;
    tdata(tindexer.translate("control1", "control2")) = 1.;
    tdata(tindexer.translate("control1", "control3")) = 1.;
    tdata(tindexer.translate("control1", "control4")) = 1.;
    tdata(tindexer.translate("control5", "control6")) = 1.;
    tdata(tindexer.translate("case1", "control2")) = 1.;
    tdata(tindexer.translate("case1", "control3")) = 1.;
    tdata(tindexer.translate("case1", "control4")) = 1.;
    tdata(tindexer.translate("case5", "control6")) = 1.;

    // Follow normal execution
    double statistic;

    int64_t cscs = 0;
    int64_t cscn = 0;
    int64_t cncn = 0;

    std::vector<std::pair<size_t, size_t>> tpairs;

    if (tpairs.empty()) {
        for (auto it = tdata.begin(); it != tdata.end(); ++it) {
            auto p = tindexer.back_translate(it.row());
            try {
                tpairs.emplace_back(p);
            } catch (std::length_error &e) {
                std::cerr << "Failed to emplace or push at " << __LINE__ << std::endl;
                throw(e);
            }
        }
    }

    for (auto &p : tpairs) {
        auto &[p1, p2] = p;
        int x = tphenotypes_[p1];
        int y = tphenotypes_[p2];

        cscs += ((x == 1) && (y == 1));
        cscn += ((x == 1) && (y == 0));
        cscn += ((x == 0) && (y == 1));
        cncn += ((x == 0) && (y == 0));
    }

    statistic = static_cast<double>(cscs) / tindexer.case_case - static_cast<double>(cscn) / tindexer.case_cont - static_cast<double>(cncn) / tindexer.cont_cont;
    fmt::print(std::cerr, "Test statistic: {}\n", statistic);
    fmt::print(std::cerr, "cscs: {}, cscn: {}, cncn: {}\n", cscs, cscn, cncn);
    fmt::print(std::cerr,
               "cscs_count: {}, cscn_count: {}, cncn_count: {}\n",
               tindexer.case_case,
               tindexer.case_cont,
               tindexer.cont_cont);
}

template class Statistic<pheno_vector>;
template class Statistic<compressed_pheno_vector>;
