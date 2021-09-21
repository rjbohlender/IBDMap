//
// Created by Bohlender,Ryan James on 11/12/20.
//

#include "source.hpp"
#include "isgzipped.hpp"
#include <boost/iostreams/filter/gzip.hpp>

/**
 * @brief Extend the lifetime of the filtering stream.
 * @param path Path to a file that may be gzipped.
 */
Source::Source(const std::string &path) {
    streambuf = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    if (is_gzipped(path)) {
        ifs.open(path, std::ios_base::in | std::ios_base::binary);
        (*streambuf).push(boost::iostreams::gzip_decompressor());
        (*streambuf).push(ifs);
    } else {
        ifs.open(path, std::ios_base::in);
        (*streambuf).push(ifs);
    }
}
