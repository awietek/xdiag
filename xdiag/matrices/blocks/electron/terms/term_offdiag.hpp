// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::electron {

// General off-diagonal term on the symmetric electron basis: the fully generic
// kernel where a single operator may move BOTH sectors at once.
// `apply(ups, dns)` returns {coeff, ups_out, dns_out}; `coeff` must already
// carry every operator-intrinsic sign (Jordan-Wigner / cross sign), as the
// kernel only adds the symmetrisation weights (character, norms, and the fermi
// sign of mapping the output product state to its representative).
//
// term_ups / term_dns are sector-restricted optimisations of this; term_diag
// forwards here for its cross-block (different irrep) case. The per-element
// representative lookup (representative_data_fermi) makes this the slow,
// general fallback, used where the sector-restricted kernels do not apply.
template <typename coeff_t, typename basis_t, typename apply_f,
          typename fill_f,
          typename = decltype(std::declval<basis_t>().dns_for_ups_rep(0))>
void term_offdiag(basis_t const &basis_in, basis_t const &basis_out,
                  apply_f apply, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  auto const &basis_up = basis_in.basis_up();
  auto bloch = basis_out.characters().template as<arma::Col<coeff_t>>();

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
      int64_t off_in = basis_in.ups_offset(idx_ups);
      auto norms_in = basis_in.norms_for_ups_rep(idx_ups);
      int64_t idx_dn = 0;
      for (bit_t dns : basis_in.dns_for_ups_rep(idx_ups)) {
        auto [c, ups_out, dns_out] = apply(ups, dns);
        auto [idx_out, sym, norm_out, fermi] =
            basis_out.representative_data_fermi(ups_out, dns_out);
        if (idx_out >= 0) {
          coeff_t val = c * bloch(sym) * norm_out / norms_in[idx_dn];
          XDIAG_FILL(off_in + idx_dn, idx_out, fermi ? -val : val);
        }
        ++idx_dn;
      }
    }
#ifdef _OPENMP
  }
#endif
}

} // namespace xdiag::matrices::electron
