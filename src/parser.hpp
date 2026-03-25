//
// Created by Bohlender,Ryan James on 9/4/19.
//

#ifndef CARVAIBD_PARSER_HPP
#define CARVAIBD_PARSER_HPP

#include <armadillo>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../binary_format/ibd_reader.hpp"
#include "breakpoint.hpp"
#include "geneticmap.hpp"
#include "indexer.hpp"
#include "math.hpp"
#include "parameters.hpp"
#include "phenotypes.hpp"
#include "reporter.hpp"
#include "statistic.hpp"
#include "threadpool.hpp"

struct Chunk {
    uint64_t start;  // index of checkpoint position
    uint64_t end;    // one past last position in this chunk
};

template <typename T>
class Parser {
    void parse_binary(IbdfReader &reader, const std::vector<std::string> &sample_names);

    std::vector<Chunk> find_chunks(const IbdfReader &reader);

    void parse_chunk(const IbdfReader &reader,
                     const std::vector<std::string> &sample_names,
                     const Chunk &chunk,
                     std::shared_ptr<Reporter> chunk_reporter,
                     ThreadPool<Statistic<T>> &threadpool);

    void concatenate_chunk_files(const std::vector<std::string> &chunk_paths);

    bool check_range(int pos);
    bool check_exclude(int pos);
    bool check_r2(const arma::SpCol<int32_t> &data, const arma::SpCol<int32_t> &last);

    void apply_checkpoint_to_sparse(arma::SpCol<int32_t> &data,
                                    const std::vector<ibdf::Segment> &active,
                                    const std::vector<std::string> &sample_names);

    void apply_delta_to_sparse(arma::SpCol<int32_t> &data,
                               const ibdf::DeltaView &dv,
                               const std::vector<std::string> &sample_names,
                               int value);

    std::vector<std::string> load_sample_names(const std::string &samples_path);

public:
    Parameters params;
    GeneticMap gmap;
    Phenotypes<T> pheno;

    Parser(Parameters params_, GeneticMap &gmap_, Phenotypes<T> &pheno_);
};

#endif//CARVAIBD_PARSER_HPP
