//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "source.hpp"
#include "iscompressed.hpp"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zstd.hpp>

/**
 * @brief Extend the lifetime of the filtering stream.
 * @param path Path to a file that may be gzipped.
 */
Source::Source(const std::string &path) {
    streambuf = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    switch (is_compressed(path)) {
        case CompressionType::gzip:
            ifs.open(path, std::ios_base::in | std::ios_base::binary);
            (*streambuf).push(boost::iostreams::gzip_decompressor());
            (*streambuf).push(ifs);
            break;
        case CompressionType::zstd:
            ifs.open(path, std::ios_base::in | std::ios_base::binary);
            (*streambuf).push(boost::iostreams::zstd_decompressor());
            (*streambuf).push(ifs);
            break;
        case CompressionType::uncompressed:
            ifs.open(path, std::ios_base::in);
            (*streambuf).push(ifs);
            break;
        default: break;
    }
}
