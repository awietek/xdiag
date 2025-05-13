// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <tuple>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_szsz(Coupling const &cpl, Op const &op, basis_t const &basis,
                arma::Col<coeff_t> const &vec_in,
                arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t ss1 = op[0];
  int64_t ss2 = op[1];
  coeff_t val_same = J / 4.;
  coeff_t val_diff = -J / 4.;

  int64_t s1 = std::min(ss1, ss2);
  int64_t s2 = std::max(ss1, ss2);

  int64_t n_postfix_bits = basis.n_postfix_bits();
  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << (s2 - n_postfix_bits);

  int64_t idx = 0;
  for (bit_t prefix : basis.prefixes()) {
    bit_t prefix_shifted = (prefix << n_postfix_bits);
    auto const &postfixes = basis.postfix_states(prefix);

    // Both sites are on prefixes
    if ((s1 >= n_postfix_bits) && (s2 >= n_postfix_bits)) {
      coeff_t val =
          (bits::popcnt(prefix_shifted & mask) & 1) ? val_diff : val_same;
      int64_t end = idx + postfixes.size();
      for (; idx < end; ++idx) {
        vec_out(idx) += val * vec_in(idx);
      }
    }

    // Both sites are on postfixes
    else if ((s1 < n_postfix_bits) && (s2 < n_postfix_bits)) {

      for (auto postfix : postfixes) {
        if (bits::popcnt(postfix & mask) & 1) {
          vec_out(idx) += val_diff * vec_in(idx);
        } else {
          vec_out(idx) += val_same * vec_in(idx);
        }
        ++idx;
      }
    }

    // s2 is prefix s1 is postfix
    else {

      // s2 is up
      if (prefix & s2mask) {
        for (auto postfix : postfixes) {
          if (postfix & s1mask) {
            vec_out(idx) += val_same * vec_in(idx);
          } else {
            vec_out(idx) += val_diff * vec_in(idx);
          }
          ++idx;
        }
      }

      // s2 is dn
      else {
        for (auto postfix : postfixes) {
          if (postfix & s1mask) {
            vec_out(idx) += val_diff * vec_in(idx);
          } else {
            vec_out(idx) += val_same * vec_in(idx);
          }
          ++idx;
        }
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf_distributed
