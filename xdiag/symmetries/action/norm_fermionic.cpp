// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "norm_fermionic.hpp"

#include <cmath>
#include <cstdint>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/action/site_permutation_sublattice.hpp>
#include <xdiag/symmetries/fermi_sign.hpp>

namespace xdiag::symmetries {

template <typename bit_t, typename coeff_t, typename action_t>
double norm_fermionic(bit_t state, action_t const &action,
                      arma::Col<coeff_t> const &characters) {
  auto const &group = action.group();
  coeff_t amplitude = 0.0;
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    if (action.apply(sym, state) == state) {
      bool fermi = fermi_bool_of_permutation(state, group[sym]);
      amplitude += fermi ? -characters(sym) : characters(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

#define INSTANTIATE_NORM_FERMIONIC_SP(BIT_TYPE)                                \
  template double norm_fermionic(BIT_TYPE, SitePermutation const &,            \
                                 arma::vec const &);                           \
  template double norm_fermionic(BIT_TYPE, SitePermutation const &,            \
                                 arma::cx_vec const &);                        \
  using namespace bits;

INSTANTIATE_NORM_FERMIONIC_SP(uint32_t);
INSTANTIATE_NORM_FERMIONIC_SP(uint64_t);
INSTANTIATE_NORM_FERMIONIC_SP(BitsetDynamic);
INSTANTIATE_NORM_FERMIONIC_SP(BitsetStatic1);
INSTANTIATE_NORM_FERMIONIC_SP(BitsetStatic2);
INSTANTIATE_NORM_FERMIONIC_SP(BitsetStatic4);
INSTANTIATE_NORM_FERMIONIC_SP(BitsetStatic8);

#undef INSTANTIATE_NORM_FERMIONIC_SP
} // namespace xdiag::symmetries
