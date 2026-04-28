#ifndef CARVAIBD_PARQUET_REPORTER_HPP
#define CARVAIBD_PARQUET_REPORTER_HPP

#include "resultrow.hpp"
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if defined(IBDMAP_HAS_PARQUET)
#include <arrow/api.h>
#include <parquet/arrow/writer.h>
#endif

class ParquetReporter {
    std::atomic_bool done;
    std::atomic_bool failed;
    std::atomic<int> nrows;
    std::atomic<int> nsubmitted;
    std::atomic<int> nwritten;
    std::string out_path;
    std::thread print_thread;
    std::map<uint64_t, ResultRow> reorder_buf;
    uint64_t next_seq = 0;
    mutable std::mutex mut;
    std::condition_variable data_cond;
    size_t batch_size;

#if defined(IBDMAP_HAS_PARQUET)
    std::shared_ptr<arrow::Schema> schema;
    std::unique_ptr<parquet::arrow::FileWriter> writer;

    arrow::StringBuilder chrom_builder;
    arrow::Int32Builder pos_builder;
    arrow::DoubleBuilder orig_cscs_rate_builder;
    arrow::DoubleBuilder orig_cscn_rate_builder;
    arrow::DoubleBuilder orig_cncn_rate_builder;
    arrow::DoubleBuilder original_builder;
    std::unique_ptr<arrow::ListBuilder> permutations_builder;
    arrow::FloatBuilder *permutation_value_builder = nullptr;
    int64_t buffered_rows = 0;

    void open_writer();
    void append_row(const ResultRow &row);
    void flush_batch();
#endif

    void print();

public:
    explicit ParquetReporter(std::string output, size_t batch_size_ = 4096);
    ~ParquetReporter();

    void submit(ResultRow row);

    static void merge_files(const std::vector<std::string> &input_paths,
                            const std::string &output_path);
};

#endif//CARVAIBD_PARQUET_REPORTER_HPP
