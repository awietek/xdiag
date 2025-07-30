// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_szsz(Coupling const &cpl, Op const &op, basis_t const &basis,
                fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  // Set values for same/diff (tJ block definition)
  coeff_t val_same = J / 4.;
  coeff_t val_diff = -J / 4.;

  // bitmasks for fast evaluations
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << s2;
  bit_t mask = s1mask | s2mask;

  if constexpr (symmetric) {

#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
      for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
        bit_t ups = basis.rep_ups(idx_ups);
        auto dnss = basis.dns_for_ups_rep(ups);
        int64_t idx = basis.ups_offset(idx_ups);

        if ((ups & mask) == mask) { // both spins pointing up
          for (bit_t dns : dnss) {
            if (!(dns & mask)) {
              XDIAG_FILL(idx, idx, val_same);
            }
            ++idx;
          }
        } else if (ups & s1mask) { // s1 is pointing up
          for (bit_t dns : dnss) {
            if ((dns & mask) == s2mask) {
              XDIAG_FILL(idx, idx, val_diff);
            }
            ++idx;
          }
        } else if (ups & s2mask) { // s2 is pointing up
          for (bit_t dns : dnss) {
            if ((dns & mask) == s1mask) {
              XDIAG_FILL(idx, idx, val_diff);
            }
            ++idx;
          }
        } else { // no upspins
          for (bit_t dns : dnss) {
            if ((dns & mask) == mask) {
              XDIAG_FILL(idx, idx, val_same);
            }
            ++idx;
          }
        }
      }
#ifdef _OPENMP
    }
#endif
  } else { // if not (symmetric)

#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
      auto ups_and_idces = basis.states_indices_ups_thread();
#else
    auto ups_and_idces = basis.states_indices_ups();
#endif

      for (auto [up, idx_up] : ups_and_idces) {

        int64_t idx = idx_up * basis.size_dns();

        if ((up & mask) == mask) { // both spins pointing up
          for (bit_t dn : basis.states_dns()) {
            if (!(dn & mask)) {
              XDIAG_FILL(idx, idx, val_same);
            }
            ++idx;
          }
        } else if (up & s1mask) { // s1 is pointing up
          for (bit_t dn : basis.states_dns()) {
            if ((dn & mask) == s2mask) {
              XDIAG_FILL(idx, idx, val_diff);
            }
            ++idx;
          }
        } else if (up & s2mask) { // s2 is pointing up
          for (bit_t dn : basis.states_dns()) {
            if ((dn & mask) == s1mask) {
              XDIAG_FILL(idx, idx, val_diff);
            }
            ++idx;
          }

        } else { // no upspins
          for (bit_t dn : basis.states_dns()) {
            if ((dn & mask) == mask) {
              XDIAG_FILL(idx, idx, val_same);
            }
            ++idx;
          }
        }

      } // for (auto up : ...)

#ifdef _OPENMP
    }
#endif
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
