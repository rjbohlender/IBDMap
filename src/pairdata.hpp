#ifndef IBDMAP_PAIRDATA_HPP
#define IBDMAP_PAIRDATA_HPP

#include <armadillo>
#include <cstdint>
#include <type_traits>

// Armadillo names its supported signed byte type arma::s8. On platforms where
// int8_t is signed char and arma::s8 is char, they are distinct C++ types even
// though both occupy one signed byte.
using pair_count_t = arma::s8;
using pair_data_t = arma::SpCol<pair_count_t>;

static_assert(sizeof(pair_count_t) == sizeof(int8_t));
static_assert(std::is_signed_v<pair_count_t>);

#endif//IBDMAP_PAIRDATA_HPP
