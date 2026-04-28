#include "parquet_reporter.hpp"
#include <chrono>
#include <fmt/include/fmt/ostream.h>
#include <iostream>
#include <stdexcept>
#include <utility>

#if defined(IBDMAP_HAS_PARQUET)
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>
#endif

#if defined(IBDMAP_HAS_PARQUET)
namespace {
void check_status(const arrow::Status &status, const std::string &context) {
    if (!status.ok()) {
        throw std::runtime_error(context + ": " + status.ToString());
    }
}

template<typename T>
T value_or_throw(arrow::Result<T> result, const std::string &context) {
    if (!result.ok()) {
        throw std::runtime_error(context + ": " + result.status().ToString());
    }
    return std::move(*result);
}

std::shared_ptr<arrow::Schema> result_schema() {
    return arrow::schema({
            arrow::field("chrom", arrow::utf8(), false),
            arrow::field("pos", arrow::int32(), false),
            arrow::field("orig_cscs_rate", arrow::float64(), false),
            arrow::field("orig_cscn_rate", arrow::float64(), false),
            arrow::field("orig_cncn_rate", arrow::float64(), false),
            arrow::field("original", arrow::float64(), false),
            arrow::field("permutations", arrow::list(arrow::float32()), false),
    });
}
}// namespace
#endif

ParquetReporter::ParquetReporter(std::string output, size_t batch_size_) :
        done(false),
        failed(false),
        nrows(0),
        nsubmitted(0),
        nwritten(0),
        out_path(std::move(output)),
        batch_size(batch_size_) {
#if defined(IBDMAP_HAS_PARQUET)
    schema = result_schema();
    permutations_builder = std::make_unique<arrow::ListBuilder>(
            arrow::default_memory_pool(),
            std::make_shared<arrow::FloatBuilder>());
    permutation_value_builder = static_cast<arrow::FloatBuilder *>(permutations_builder->value_builder());
    print_thread = std::thread(&ParquetReporter::print, std::ref(*this));
#else
    throw std::runtime_error("Parquet output was requested, but this build does not include Arrow/Parquet support.");
#endif
}

ParquetReporter::~ParquetReporter() {
    while (!failed && (!reorder_buf.empty() || nrows > 0)) {
        fmt::print(std::cerr, "parquet queue size: {} nrows: {}\n", reorder_buf.size(), nrows);
        std::this_thread::sleep_for(std::chrono::nanoseconds(100000000));
    }
    done = true;
    data_cond.notify_all();
    if (print_thread.joinable()) {
        print_thread.join();
    }
}

void ParquetReporter::submit(ResultRow row) {
    if (failed) {
        throw std::runtime_error("cannot submit to failed ParquetReporter");
    }
    std::unique_lock<std::mutex> lk(mut);
    reorder_buf[row.seq] = std::move(row);
    nrows++;
    nsubmitted++;
    lk.unlock();
    data_cond.notify_all();
}

void ParquetReporter::print() {
#if defined(IBDMAP_HAS_PARQUET)
    try {
        open_writer();
        std::unique_lock<std::mutex> lk(mut);
        while (!done || nrows > 0) {
            data_cond.wait_for(lk, std::chrono::seconds(1), [this] {
                return reorder_buf.count(next_seq) || done;
            });
            while (reorder_buf.count(next_seq)) {
                auto node = reorder_buf.extract(next_seq);
                lk.unlock();
                append_row(node.mapped());
                if (buffered_rows >= static_cast<int64_t>(batch_size)) {
                    flush_batch();
                }
                nrows--;
                nwritten++;
                lk.lock();
                next_seq++;
            }
        }
        lk.unlock();
        flush_batch();
        check_status(writer->Close(), "closing parquet writer");
    } catch (const std::exception &e) {
        fmt::print(std::cerr, "ParquetReporter failed: {}\n", e.what());
        failed = true;
        nrows = 0;
        done = true;
    }
#endif
}

#if defined(IBDMAP_HAS_PARQUET)
void ParquetReporter::open_writer() {
    auto outfile = value_or_throw(arrow::io::FileOutputStream::Open(out_path), "opening parquet output");
    auto props = parquet::WriterProperties::Builder()
            .compression(parquet::Compression::ZSTD)
            ->build();
    writer = value_or_throw(
            parquet::arrow::FileWriter::Open(
                    *schema,
                    arrow::default_memory_pool(),
                    outfile,
                    props,
                    parquet::ArrowWriterProperties::Builder().build()),
            "opening parquet writer");
}

void ParquetReporter::append_row(const ResultRow &row) {
    check_status(chrom_builder.Append(row.chrom), "appending chrom");
    check_status(pos_builder.Append(row.pos), "appending pos");
    check_status(orig_cscs_rate_builder.Append(row.orig_cscs_rate), "appending orig_cscs_rate");
    check_status(orig_cscn_rate_builder.Append(row.orig_cscn_rate), "appending orig_cscn_rate");
    check_status(orig_cncn_rate_builder.Append(row.orig_cncn_rate), "appending orig_cncn_rate");
    check_status(original_builder.Append(row.original), "appending original");
    check_status(permutations_builder->Append(), "appending permutations list");
    check_status(permutation_value_builder->AppendValues(row.permutations), "appending permutation values");
    buffered_rows++;
}

void ParquetReporter::flush_batch() {
    if (buffered_rows == 0) {
        return;
    }

    std::shared_ptr<arrow::Array> chrom;
    std::shared_ptr<arrow::Array> pos;
    std::shared_ptr<arrow::Array> orig_cscs_rate;
    std::shared_ptr<arrow::Array> orig_cscn_rate;
    std::shared_ptr<arrow::Array> orig_cncn_rate;
    std::shared_ptr<arrow::Array> original;
    std::shared_ptr<arrow::Array> permutations;

    check_status(chrom_builder.Finish(&chrom), "finishing chrom array");
    check_status(pos_builder.Finish(&pos), "finishing pos array");
    check_status(orig_cscs_rate_builder.Finish(&orig_cscs_rate), "finishing orig_cscs_rate array");
    check_status(orig_cscn_rate_builder.Finish(&orig_cscn_rate), "finishing orig_cscn_rate array");
    check_status(orig_cncn_rate_builder.Finish(&orig_cncn_rate), "finishing orig_cncn_rate array");
    check_status(original_builder.Finish(&original), "finishing original array");
    check_status(permutations_builder->Finish(&permutations), "finishing permutations array");

    auto batch = arrow::RecordBatch::Make(
            schema,
            buffered_rows,
            {chrom, pos, orig_cscs_rate, orig_cscn_rate, orig_cncn_rate, original, permutations});
    check_status(writer->WriteRecordBatch(*batch), "writing parquet record batch");

    chrom_builder.Reset();
    pos_builder.Reset();
    orig_cscs_rate_builder.Reset();
    orig_cscn_rate_builder.Reset();
    orig_cncn_rate_builder.Reset();
    original_builder.Reset();
    permutations_builder = std::make_unique<arrow::ListBuilder>(
            arrow::default_memory_pool(),
            std::make_shared<arrow::FloatBuilder>());
    permutation_value_builder = static_cast<arrow::FloatBuilder *>(permutations_builder->value_builder());
    buffered_rows = 0;
}
#endif

void ParquetReporter::merge_files(const std::vector<std::string> &input_paths,
                                  const std::string &output_path) {
#if defined(IBDMAP_HAS_PARQUET)
    auto schema = result_schema();
    auto outfile = value_or_throw(arrow::io::FileOutputStream::Open(output_path), "opening merged parquet output");
    auto props = parquet::WriterProperties::Builder()
            .compression(parquet::Compression::ZSTD)
            ->build();
    auto writer = value_or_throw(
            parquet::arrow::FileWriter::Open(
                    *schema,
                    arrow::default_memory_pool(),
                    outfile,
                    props,
                    parquet::ArrowWriterProperties::Builder().build()),
            "opening merged parquet writer");

    for (const auto &path : input_paths) {
        auto infile = value_or_throw(arrow::io::ReadableFile::Open(path), "opening parquet chunk");
        auto reader = value_or_throw(
                parquet::arrow::OpenFile(infile, arrow::default_memory_pool()),
                "opening parquet reader");

        auto batch_reader = value_or_throw(
                reader->GetRecordBatchReader(),
                "creating parquet batch reader");

        while (true) {
            auto maybe_batch = batch_reader->Next();
            if (!maybe_batch.ok()) {
                throw std::runtime_error("reading parquet record batch: " + maybe_batch.status().ToString());
            }
            std::shared_ptr<arrow::RecordBatch> batch = *maybe_batch;
            if (!batch) {
                break;
            }
            check_status(writer->WriteRecordBatch(*batch), "writing merged parquet batch");
        }
    }

    check_status(writer->Close(), "closing merged parquet writer");
#else
    (void) input_paths;
    (void) output_path;
    throw std::runtime_error("Parquet merge requested, but this build does not include Arrow/Parquet support.");
#endif
}
