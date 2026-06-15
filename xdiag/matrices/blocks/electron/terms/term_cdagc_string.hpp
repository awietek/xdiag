// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/matrices/terms/cdagc_string.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::electron {

// Applies a product of elementary electron operators (Cdagup/Cup/Cdagdn/Cdn) to
// the electron product basis -- the spinful analogue of
// fermion::term_cdagc_string. The electron block is two independent fermion
// sectors in the "all ups then all dns" Jordan-Wigner order, so the string
// splits into an up sub-string acting on ups and a dn sub-string acting on dns,
// each handled by the shared CdagCString. The only spinful addition is the
// cross sign: every dn operator sits behind the whole up string, so the dn
// block carries (-1)^(q_dn * Nup) with q_dn the number of dn operators and Nup
// the up occupation (evaluated on the incoming ups, since in this order the dn
// operators act first).
//
// The up sub-string and the per-up cross sign are evaluated once per up state in
// the outer loop; the dn sub-string action is precomputed once per dn state.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagc_string(Coeff const &c, Monomial const &mono,
                       basis::BasisElectron<enumeration_t> const &basis_in,
                       basis::BasisElectron<enumeration_t> const &basis_out,
                       fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();

  // Split the monomial into its up and dn sub-strings (validating the types).
  std::vector<Op> up_ops, dn_ops;
  for (int64_t k = 0; k < mono.size(); ++k) {
    Op const &op = mono[k];
    std::string type = op.type();
    if ((type == "Cdagup") || (type == "Cup")) {
      up_ops.push_back(op);
    } else if ((type == "Cdagdn") || (type == "Cdn")) {
      dn_ops.push_back(op);
    } else {
      XDIAG_THROW("electron Cdag/C string: unexpected operator type \"" + type +
                  "\"; only Cdagup, Cup, Cdagdn, Cdn are allowed");
    }
  }
  CdagCString<bit_t> up_str(nsites, Monomial(up_ops), "Cdagup", "Cup");
  CdagCString<bit_t> dn_str(nsites, Monomial(dn_ops), "Cdagdn", "Cdn");
  bool cross_when_odd_nup = (dn_ops.size() & 1);

  auto const &basis_up_in = basis_in.basis_up();
  auto const &basis_up_out = basis_out.basis_up();
  auto const &basis_dn_in = basis_in.basis_dn();
  auto const &basis_dn_out = basis_out.basis_dn();
  int64_t size_dn_in = basis_dn_in.size();
  int64_t size_dn_out = basis_dn_out.size();

  // Precompute the dn sub-string action for every dn state (independent of the
  // up sector). dn_out[idx] == -1 marks a dn state annihilated by the string.
  // std::vector<char> (not <bool>) avoids a data race on the packed bits.
  std::vector<int64_t> dn_out(size_dn_in);
  std::vector<char> dn_neg(size_dn_in, 0);
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_dn, end_dn, idx_dn] =
        utils::thread_range(basis_dn_in, num_thread, omp_get_num_threads());
    for (auto it = begin_dn; it != end_dn; ++it, ++idx_dn) {
      bit_t dns = *it;
      if (dn_str.non_zero(dns)) {
        std::pair<bit_t, bool> a = dn_str.action(dns);
        dn_out[idx_dn] = basis_dn_out.index(a.first);
        dn_neg[idx_dn] = a.second ? 1 : 0;
      } else {
        dn_out[idx_dn] = -1;
      }
    }
  }
#else
  {
    int64_t idx_dn = 0;
    for (bit_t dns : basis_dn_in) {
      if (dn_str.non_zero(dns)) {
        std::pair<bit_t, bool> a = dn_str.action(dns);
        dn_out[idx_dn] = basis_dn_out.index(a.first);
        dn_neg[idx_dn] = a.second ? 1 : 0;
      } else {
        dn_out[idx_dn] = -1;
      }
      ++idx_dn;
    }
  }
#endif

  // Apply the up sub-string once per up state and fold in the dn table.
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up_in, num_thread, omp_get_num_threads());
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      if (!up_str.non_zero(ups)) {
        continue;
      }
      std::pair<bit_t, bool> a = up_str.action(ups);
      bool up_total = a.second;
      if (cross_when_odd_nup && (bits::popcount(ups) & 1)) {
        up_total = !up_total;
      }
      int64_t base_in = idx_up * size_dn_in;
      int64_t base_out = basis_up_out.index(a.first) * size_dn_out;
      for (int64_t idx_dn = 0; idx_dn < size_dn_in; ++idx_dn) {
        if (dn_out[idx_dn] >= 0) {
          bool neg = up_total ^ (bool)dn_neg[idx_dn];
          XDIAG_FILL(base_in + idx_dn, base_out + dn_out[idx_dn],
                     neg ? -cf : cf);
        }
      }
    }
  }
#else
  int64_t idx_up = 0;
  for (bit_t ups : basis_up_in) {
    if (up_str.non_zero(ups)) {
      std::pair<bit_t, bool> a = up_str.action(ups);
      bool up_total = a.second;
      if (cross_when_odd_nup && (bits::popcount(ups) & 1)) {
        up_total = !up_total;
      }
      int64_t base_in = idx_up * size_dn_in;
      int64_t base_out = basis_up_out.index(a.first) * size_dn_out;
      for (int64_t idx_dn = 0; idx_dn < size_dn_in; ++idx_dn) {
        if (dn_out[idx_dn] >= 0) {
          bool neg = up_total ^ (bool)dn_neg[idx_dn];
          fill(base_in + idx_dn, base_out + dn_out[idx_dn], neg ? -cf : cf);
        }
      }
    }
    ++idx_up;
  }
#endif
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
