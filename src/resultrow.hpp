#ifndef CARVAIBD_RESULTROW_HPP
#define CARVAIBD_RESULTROW_HPP

#include <cstdint>
#include <string>
#include <vector>

struct ResultRow {
    uint64_t seq = 0;
    std::string chrom;
    int32_t pos = 0;
    double orig_cscs_rate = 0;
    double orig_cscn_rate = 0;
    double orig_cncn_rate = 0;
    double original = 0;
    std::vector<float> permutations;
};

#endif//CARVAIBD_RESULTROW_HPP
