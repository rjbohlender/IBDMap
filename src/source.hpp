//
// Created by Bohlender,Ryan James on 11/12/20.
//

#ifndef CARVAIBD_SOURCE_HPP
#define CARVAIBD_SOURCE_HPP

#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <fstream>
#include <memory>

struct Source {
    std::ifstream ifs;
    std::unique_ptr<boost::iostreams::filtering_streambuf<boost::iostreams::input>> streambuf;

    explicit Source(const std::string &path);
};

#endif//CARVAIBD_SOURCE_HPP
