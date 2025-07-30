// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_number(Coupling const &cpl, Op const &op, basis_t const &basis,
                  fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;
  std::string type = op.type();

  if (type == "Nup") {

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
          int64_t size_dns = dnss.size();
          if (ups & mask) { // check whether bit at s is set
            int64_t end = idx + size_dns;
            for (; idx < end; ++idx) {
              XDIAG_FILL(idx, idx, mu);
            }
          } else {
            idx += size_dns;
          }
        }
#ifdef _OPENMP
      }
#endif

    } else { // if not symmetric

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

          int64_t size_dns = basis.size_dns();
          if (up & mask) { // check whether bit at s is set
            int64_t end = idx + size_dns;
            for (; idx < end; ++idx) {
              XDIAG_FILL(idx, idx, mu);
            }
          } else {
            idx += size_dns;
          }
        }
#ifdef _OPENMP
      }
#endif
    }
  } else if (type == "Ndn") {

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
          for (bit_t dns : dnss) {
            if (dns & mask) {
              XDIAG_FILL(idx, idx, mu);
            }
            ++idx;
          }
        }
#ifdef _OPENMP
      }
#endif

    } else { // if not symmetric
#ifdef _OPENMP
#pragma omp parallel
      {
        int num_thread = omp_get_thread_num();
        auto ups_and_idces = basis.states_indices_ups_thread();
#else
      auto ups_and_idces = basis.states_indices_ups();
#endif

        for (auto [up, idx_up] : ups_and_idces) {
          (void)up;
          int64_t idx = idx_up * basis.size_dns();

          for (bit_t dn : basis.states_dns()) {
            if (dn & mask) {
              XDIAG_FILL(idx, idx, mu);
            }
            ++idx;
          }
        }

#ifdef _OPENMP
      }
#endif
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
