//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "parser.hpp"
#include <atomic>
#include <cstdio>
#include <fstream>
#include <random>
#include <thread>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <fmt/include/fmt/ostream.h>

template <typename T>
std::vector<std::string> Parser<T>::load_sample_names(const std::string &samples_path) {
    std::vector<std::string> names;
    std::ifstream ifs(samples_path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open samples file: " + samples_path);
    }
    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty()) {
            names.push_back(line);
        }
    }
    return names;
}

template <typename T>
void Parser<T>::apply_checkpoint_to_sparse(arma::SpCol<int32_t> &data,
                                           const std::vector<ibdf::Segment> &active,
                                           const std::vector<std::string> &sample_names) {
    data.zeros();
    for (const auto &seg : active) {
        if (params.cM > 0 && seg.cm < params.cM) {
            continue;
        }
        const auto &name1 = sample_names[seg.p1];
        const auto &name2 = sample_names[seg.p2];
        arma::sword row_idx = (*pheno.indexer).translate(name1, name2);
        if (row_idx >= 0) {
            data(row_idx) = 1;
        }
    }
}

template <typename T>
void Parser<T>::apply_delta_to_sparse(arma::SpCol<int32_t> &data,
                                      const ibdf::DeltaView &dv,
                                      const std::vector<std::string> &sample_names,
                                      int value) {
    uint32_t n = (value > 0) ? dv.n_adds() : dv.n_dels();
    const auto cm_arr  = (value > 0) ? dv.add_cm : dv.del_cm;
    const auto p1_arr  = (value > 0) ? dv.add_p1 : dv.del_p1;
    const auto p2_arr  = (value > 0) ? dv.add_p2 : dv.del_p2;

    for (uint32_t i = 0; i < n; ++i) {
        if (params.cM > 0 && cm_arr[i] < params.cM) {
            continue;
        }
        const auto &name1 = sample_names[p1_arr[i]];
        const auto &name2 = sample_names[p2_arr[i]];
        arma::sword row_idx = (*pheno.indexer).translate(name1, name2);
        if (row_idx >= 0) {
            data(row_idx) += value;
        }
    }
}

template <typename T>
std::vector<Chunk> Parser<T>::find_chunks(const IbdfReader &reader) {
    std::vector<Chunk> chunks;
    for (uint64_t i = 0; i < reader.n_positions(); ++i) {
        if (reader.index()[i].is_checkpoint()) {
            if (!chunks.empty()) {
                chunks.back().end = i;
            }
            chunks.push_back({i, reader.n_positions()});
        }
    }
    return chunks;
}

template <typename T>
void Parser<T>::parse_chunk(const IbdfReader &reader,
                            const std::vector<std::string> &sample_names,
                            const Chunk &chunk,
                            std::shared_ptr<Reporter> chunk_reporter,
                            ThreadPool<Statistic<T>> &threadpool) {
    arma::uword n_pairs = pheno.samples->size() * (pheno.samples->size() - 1) / 2;
    arma::SpCol<int32_t> data(n_pairs);
    arma::SpCol<int32_t> last(n_pairs);

    double cur_dist = 0;
    double last_dist = 0;
    uint64_t seq = 0;

    std::vector<ibdf::Segment> active;

    for (uint64_t i = chunk.start; i < chunk.end; ++i) {
        const auto &idx_entry = reader.index()[i];
        uint64_t bp_pos = idx_entry.bp_pos;
        int pos = static_cast<int>(bp_pos);

        const auto *hdr = ibdf::peek_header(reader.block_ptr(i));
        if (hdr->block_type == ibdf::BLOCK_CHECKPOINT) {
            auto cp = ibdf::decode_checkpoint(reader.block_ptr(i));
            ibdf::load_checkpoint(active, cp);
            apply_checkpoint_to_sparse(data, active, sample_names);
        } else {
            auto dv = ibdf::decode_delta(reader.block_ptr(i));
            apply_delta_to_sparse(data, dv, sample_names, -1);
            apply_delta_to_sparse(data, dv, sample_names, 1);
        }

        std::string chrom = params.chromosome;
        std::pair<std::pair<int, double>, std::pair<int, double>> nearest = gmap.find_nearest(chrom, pos);
        if (nearest.first.first != nearest.second.first) {
            cur_dist = (pos - nearest.first.first) * (nearest.second.second - nearest.first.second) / (nearest.second.first - nearest.first.first) + nearest.first.second;
        } else {
            cur_dist = nearest.first.second;
        }

        if (last_dist != 0 && cur_dist - last_dist < params.min_dist) {
            continue;
        }

        if (params.range) {
            if (check_range(pos)) {
                continue;
            }
        }
        if (params.exclude) {
            if (check_exclude(pos)) {
                fmt::print(std::cerr, "Excluded breakpoint: {}\n", pos);
                continue;
            }
        }

        if (params.rsquared) {
            if (check_r2(data, last)) {
                continue;
            }
        }
        if (arma::accu(data) < params.threshold) {
            continue;
        }

        last_dist = cur_dist;
        last = data;

        Breakpoint bp{};
        bp.breakpoint = std::make_pair(chrom, std::to_string(pos));

        Statistic stat(data,
                       bp,
                       pheno.indexer,
                       chunk_reporter,
                       seq++,
                       params,
                       pheno.phenotypes,
                       pheno.transposed);
        threadpool.submit(std::move(stat));
    }
}

template <typename T>
void Parser<T>::concatenate_chunk_files(const std::vector<std::string> &chunk_paths) {
    boost::iostreams::filtering_ostream os;
    boost::iostreams::file_sink sink(params.output_path);
    os.push(boost::iostreams::zstd_compressor());
    os.push(sink);

    for (const auto &path : chunk_paths) {
        boost::iostreams::filtering_istream is;
        boost::iostreams::file_source source(path);
        is.push(boost::iostreams::zstd_decompressor());
        is.push(source);

        std::string line;
        while (std::getline(is, line)) {
            fmt::print(os, "{}\n", line);
        }

        std::remove(path.c_str());
    }
}

template <typename T>
void Parser<T>::parse_binary(IbdfReader &reader, const std::vector<std::string> &sample_names) {
    auto chunks = find_chunks(reader);
    if (chunks.empty()) return;

    unsigned n_parse_threads = std::min(
        static_cast<unsigned>(chunks.size()),
        std::max(1u, static_cast<unsigned>(params.nthreads) / 4));

    ThreadPool<Statistic<T>> threadpool(params);

    // Generate unique temp file paths per chunk
    std::random_device rd;
    uint32_t job_id = rd();
    std::vector<std::string> chunk_paths(chunks.size());
    for (size_t i = 0; i < chunks.size(); ++i) {
        chunk_paths[i] = fmt::format("{}.chunk_{}_{:08x}.zst",
                                     params.output_path, i, job_id);
    }

    // Create reporters up front so they're alive for the duration of parsing
    std::vector<std::shared_ptr<Reporter>> chunk_reporters(chunks.size());
    for (size_t i = 0; i < chunks.size(); ++i) {
        chunk_reporters[i] = std::make_shared<Reporter>(chunk_paths[i], false);
    }

    std::atomic<size_t> next_chunk{0};
    std::vector<std::thread> parse_threads;

    for (unsigned t = 0; t < n_parse_threads; ++t) {
        parse_threads.emplace_back([&]() {
            while (true) {
                size_t ci = next_chunk.fetch_add(1);
                if (ci >= chunks.size()) break;
                parse_chunk(reader, sample_names, chunks[ci],
                            chunk_reporters[ci], threadpool);
            }
        });
    }

    for (auto &t : parse_threads) {
        t.join();
    }
    threadpool.wait_for_completion();

    // Flush and close all chunk reporters after all stats have completed
    chunk_reporters.clear();

    if (!params.output_path.empty()) {
        concatenate_chunk_files(chunk_paths);
    }
}

template <typename T>
bool Parser<T>::check_range(int pos) {
    return pos < (*params.range)[0] || pos > (*params.range)[1];
}

template <typename T>
bool Parser<T>::check_exclude(int pos) {
    bool exclude_region = false;
    for (const auto &v : *params.exclude) {
        exclude_region |= pos >= v.first && pos <= v.second;
    }
    return exclude_region;
}

template <typename T>
bool Parser<T>::check_r2(const arma::SpCol<int32_t> &data, const arma::SpCol<int32_t> &last) {
    double r2 = cor(data, last);
    if (params.verbose) {
        std::cerr << "r2: " << r2 << std::endl;
    }
    return r2 > *params.rsquared;
}

template <typename T>
Parser<T>::Parser(Parameters params_, GeneticMap &gmap_, Phenotypes<T> &pheno_)
    : params(std::move(params_)), gmap(std::move(gmap_)), pheno(pheno_) {

    std::string samples_path = params.samples;
    std::vector<std::string> sample_names = load_sample_names(samples_path);

    IbdfReader reader(params.input.c_str());

    if (params.verbose) {
        fmt::print(std::cerr, "IBDF v3: {} positions, {} samples, ckpt every {}\n",
                   reader.n_positions(), reader.n_samples(), reader.checkpoint_interval());
    }

    if (params.verbose) {
        std::cerr << "Parsing binary data\n";
    }
    parse_binary(reader, sample_names);
}

template class Parser<pheno_vector>;
template class Parser<compressed_pheno_vector>;
