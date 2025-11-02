// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t, class basis_t, class term_action_f, class fill_f>
void generic_term_diag_sym(basis_t const &basis, term_action_f term_action,
                           fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis.nsites();
  bit_t sitesmask = ((bit_t)1 << nsites) - 1;

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
    for (int64_t idx_up_in = 0; idx_up_in < basis.n_rep_ups(); ++idx_up_in) {
      int64_t idx = basis.ups_offset(idx_up_in);

      bit_t up_in = basis.rep_ups(idx_up_in);
      bit_t not_up_in = (~up_in) & sitesmask;

      auto dnss_in = basis.dns_for_ups_rep(up_in);
      auto up_syms = basis.syms_ups(up_in);

      // Trivial stabilizer of ups -> dns need to be deposited
      if (up_syms.size() == 1) {
        for (bit_t dnc_in : dnss_in) {
          int64_t dn_in = bits::deposit(dnc_in, not_up_in);
          coeff_t val = term_action(up_in, dn_in);
          XDIAG_FILL(idx, idx, val);
          ++idx;
        }
      }

      // stabilizer for ups -> dns don't need to be deposited
      else {
        for (bit_t dn_in : dnss_in) {
          coeff_t val = term_action(up_in, dn_in);
          XDIAG_FILL(idx, idx, val);
          ++idx;
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}

template <typename coeff_t, class basis_t, class term_action_f, class fill_f>
void generic_term_diag_no_sym(basis_t const &basis, term_action_f term_action,
                              fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis.nsites();
  bit_t sitesmask = ((bit_t)1 << nsites) - 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto ups_and_idces = basis.states_indices_ups_thread();
#else
  auto ups_and_idces = basis.states_indices_ups();
#endif
    for (auto [up_in, idx_up_in] : ups_and_idces) {
      bit_t not_up_in = (~up_in) & sitesmask;
      auto dncs_in = basis.states_dncs(up_in);
      int64_t idx = basis.ups_offset(idx_up_in);
      for (bit_t dnc_in : dncs_in) {
        bit_t dn_in = bits::deposit(dnc_in, not_up_in);
        coeff_t val = term_action(up_in, dn_in);
        XDIAG_FILL(idx, idx, val);
        ++idx;
      }
    }
#ifdef _OPENMP
  }
#endif
}

template <bool symmetric, typename coeff_t, class basis_t, class term_action_f,
          class fill_f>
void generic_term_diag(basis_t const &basis_in, basis_t const &basis_out,
                       term_action_f term_action, fill_f fill) {
  if (basis_in == basis_out) {
    if constexpr (symmetric) {
      generic_term_diag_sym<coeff_t>(basis_in, term_action, fill);
    } else {
      generic_term_diag_no_sym<coeff_t>(basis_in, term_action, fill);
    }
  } else {
    using bit_t = typename basis_t::bit_t;
    auto non_zero_term_ups = [&](bit_t up) { return true; };
    auto non_zero_term_dns = [&](bit_t dn) { return true; };
    auto term_actionups = [&](bit_t up) -> bit_t { return up; };
    auto term_actiondns = [&](bit_t dn) -> bit_t { return dn; };
    auto matrix_element = [&](bit_t up, bit_t dn) -> coeff_t {
      return term_action(up, dn);
    };
    generic_term_mixed<symmetric, coeff_t>(
        basis_in, basis_out, non_zero_term_ups, non_zero_term_dns,
        term_actionups, term_actiondns, matrix_element, fill);
  }
}

} // namespace xdiag::basis::tj
