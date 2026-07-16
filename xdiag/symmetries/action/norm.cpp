// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "norm.hpp"

#include <cmath>
#include <cstdint>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/action/site_permutation_sublattice.hpp>

namespace xdiag::symmetries {

template <typename bit_t, typename coeff_t, typename action_t>
double norm(bit_t state, action_t const &action,
            arma::Col<coeff_t> const &characters) {
  coeff_t amplitude = 0.0;
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    if (action.apply(sym, state) == state) {
      amplitude += characters(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

// -- SitePermutation instantiations ------------------------------------------

#define INSTANTIATE_NORM_SP(BIT_TYPE)                                          \
  template double norm(BIT_TYPE, SitePermutation const &, arma::vec const &); \
  template double norm(BIT_TYPE, SitePermutation const &,                      \
                       arma::cx_vec const &);                                  \
  using namespace bits;

INSTANTIATE_NORM_SP(uint32_t);
INSTANTIATE_NORM_SP(uint64_t);
INSTANTIATE_NORM_SP(BitsetDynamic);
INSTANTIATE_NORM_SP(BitsetStatic1);
INSTANTIATE_NORM_SP(BitsetStatic2);
INSTANTIATE_NORM_SP(BitsetStatic4);
INSTANTIATE_NORM_SP(BitsetStatic8);

INSTANTIATE_NORM_SP(BitArray1);
INSTANTIATE_NORM_SP(BitArray2);
INSTANTIATE_NORM_SP(BitArray3);
INSTANTIATE_NORM_SP(BitArray4);
INSTANTIATE_NORM_SP(BitArray8);

INSTANTIATE_NORM_SP(BitArrayLong1);
INSTANTIATE_NORM_SP(BitArrayLong2);
INSTANTIATE_NORM_SP(BitArrayLong3);
INSTANTIATE_NORM_SP(BitArrayLong4);
INSTANTIATE_NORM_SP(BitArrayLong8);

#undef INSTANTIATE_NORM_SP

// -- SitePermutationSublattice instantiations --------------------------------

#define INSTANTIATE_NORM_SPS(BIT_TYPE, N_SUBLAT)                               \
  template double norm(BIT_TYPE,                                               \
                       SitePermutationSublattice<BIT_TYPE, N_SUBLAT> const &,  \
                       arma::vec const &);                                     \
  template double norm(BIT_TYPE,                                               \
                       SitePermutationSublattice<BIT_TYPE, N_SUBLAT> const &,  \
                       arma::cx_vec const &);

INSTANTIATE_NORM_SPS(uint32_t, 1);
INSTANTIATE_NORM_SPS(uint32_t, 2);
INSTANTIATE_NORM_SPS(uint32_t, 3);
INSTANTIATE_NORM_SPS(uint32_t, 4);
INSTANTIATE_NORM_SPS(uint32_t, 5);

INSTANTIATE_NORM_SPS(uint64_t, 1);
INSTANTIATE_NORM_SPS(uint64_t, 2);
INSTANTIATE_NORM_SPS(uint64_t, 3);
INSTANTIATE_NORM_SPS(uint64_t, 4);
INSTANTIATE_NORM_SPS(uint64_t, 5);

#undef INSTANTIATE_NORM_SPS

} // namespace xdiag::symmetries
