// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>

#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::kernels::electron {

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
template <typename coeff_t, typename basis_t, typename apply_f, typename fill_f,
          typename = decltype(std::declval<basis_t>().dns_for_ups_rep(0))>
void term_offdiag(basis_t const &basis_in, basis_t const &basis_out,
                  apply_f apply, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  auto const &basis_up = basis_in.basis_up();

  // The projected states use conj(chi); the kernels run trans(ops), which gives
  // the same matrix with the characters unconjugated (see compile.hpp)
  arma::Col<coeff_t> characters =
      basis_out.characters().template as<arma::Col<coeff_t>>();

  utils::for_each_thread_range(
      basis_up, [&](auto begin_up, auto end_up, int64_t idx_ups) {
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
              coeff_t val = c * characters(sym) * norm_out / norms_in[idx_dn];
              fill(off_in + idx_dn, idx_out, fermi ? -val : val);
            }
            ++idx_dn;
          }
        }
      });
}

} // namespace xdiag::kernels::electron
