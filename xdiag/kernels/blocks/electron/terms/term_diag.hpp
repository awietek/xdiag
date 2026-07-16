// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <tuple>
#include <utility>

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/kernels/blocks/electron/terms/term_offdiag.hpp>
#include <xdiag/kernels/fill_functions.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::kernels::electron {

// Diagonal term on the (non-symmetric) electron product basis. `apply(ups, dns)`
// returns the diagonal coefficient for the product state (ups, dns). The linear
// index is idx_up * size_dn + idx_dn (up sector outer, dn sector inner),
// matching BasisElectron::index.
//
// A diagonal operator conserves both Nup and Ndn, so it only connects a block
// to itself: when basis_in != basis_out (different number sectors) there is no
// matrix element and we return immediately.
template <typename coeff_t, typename enumeration_t, typename apply_f,
          typename fill_f>
void term_diag(basis::BasisElectron<enumeration_t> const &basis_in,
               basis::BasisElectron<enumeration_t> const &basis_out,
               apply_f apply, fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;

  if (basis_in != basis_out) {
    return;
  }

  auto const &basis_up = basis_in.basis_up();
  auto const &basis_dn = basis_in.basis_dn();
  int64_t size_dn = basis_dn.size();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up, num_thread, omp_get_num_threads());
#else
    auto [begin_up, end_up, idx_up] = utils::thread_range(basis_up, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      int64_t idx = idx_up * size_dn;
      for (bit_t dns : basis_dn) {
        coeff_t c = apply(ups, dns);
        XDIAG_FILL(idx, idx, c);
        ++idx;
      }
    }
#ifdef _OPENMP
  }
#endif
}

// Symmetric overload. When basis_in == basis_out the operator is genuinely
// diagonal in the symmetry-adapted basis (a group-commuting diagonal operator
// keeps each representative state, so the matrix element is apply(ups_rep,
// dns_rep) at the state's own index, with no character/norm factors). The linear
// index runs over up representatives (outer) and their dn block (inner):
// ups_offset(idx_ups) + idx_dn.
//
// When basis_in != basis_out the operator is still diagonal in the PRODUCT
// state, but the two symmetric bases differ (same Nup/Ndn, different irrep), so
// it is off-diagonal in the index (used e.g. for cross-sector form factors).
// Unlike the non-symmetric overload -- where basis_in != basis_out can only mean
// different particle numbers (genuinely zero) -- this is not always zero, so it
// forwards to the general term_offdiag with a state-preserving action.
template <typename coeff_t, typename basis_t, typename apply_f,
          typename fill_f,
          typename = decltype(std::declval<basis_t>().dns_for_ups_rep(0))>
void term_diag(basis_t const &basis_in, basis_t const &basis_out, apply_f apply,
               fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  if (basis_in != basis_out) {
    term_offdiag<coeff_t>(
        basis_in, basis_out,
        [&](bit_t ups, bit_t dns) -> std::tuple<coeff_t, bit_t, bit_t> {
          return {apply(ups, dns), ups, dns};
        },
        fill);
    return;
  }
  auto const &basis_up = basis_in.basis_up();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_ups] =
        utils::thread_range(basis_up, num_thread, omp_get_num_threads());
#else
    auto [begin_up, end_up, idx_ups] = utils::thread_range(basis_up, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_ups) {
      bit_t ups = *it_up;
      int64_t idx = basis_in.ups_offset(idx_ups);
      for (bit_t dns : basis_in.dns_for_ups_rep(idx_ups)) {
        coeff_t c = apply(ups, dns);
        XDIAG_FILL(idx, idx, c);
        ++idx;
      }
    }
#ifdef _OPENMP
  }
#endif
}

} // namespace xdiag::kernels::electron
