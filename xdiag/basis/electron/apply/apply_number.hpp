// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/electron/apply/generic_term_diag.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>
#include <xdiag/bits/bitops.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, typename basis_t, typename fill_f>
void apply_nup_sym(Coupling const &cpl, Op const &op, basis_t const &basis,
                   fill_f fill) {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;

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
}

template <typename coeff_t, typename basis_t, typename fill_f>
void apply_nup_no_sym(Coupling const &cpl, Op const &op, basis_t const &basis,
                      fill_f fill) {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;

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

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_nup(Coupling const &cpl, Op const &op, basis_t const &basis_in,
               basis_t const &basis_out, fill_f fill) {
  if (basis_in == basis_out) {
    if constexpr (symmetric) {
      apply_nup_sym<coeff_t>(cpl, op, basis_in, fill);
    } else {
      apply_nup_no_sym<coeff_t>(cpl, op, basis_in, fill);
    }
  } else {
    using bit_t = typename basis_t::bit_t;
    coeff_t mu = cpl.scalar().as<coeff_t>();
    int64_t s = op[0];
    auto apply = [&](bit_t ups, bit_t dns) {
      return bits::gbit(ups, s) ? mu : 0.;
    };
    generic_term_diag<symmetric, coeff_t>(basis_in, basis_out, apply, fill);
  }
}

template <typename coeff_t, typename basis_t, typename fill_f>
void apply_ndn_sym(Coupling const &cpl, Op const &op, basis_t const &basis,
                   fill_f fill) {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;

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
}

template <typename coeff_t, typename basis_t, typename fill_f>
void apply_ndn_no_sym(Coupling const &cpl, Op const &op, basis_t const &basis,
                      fill_f fill) {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;

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

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_ndn(Coupling const &cpl, Op const &op, basis_t const &basis_in,
               basis_t const &basis_out, fill_f fill) {
  if (basis_in == basis_out) {
    if constexpr (symmetric) {
      apply_ndn_sym<coeff_t>(cpl, op, basis_in, fill);
    } else {
      apply_ndn_no_sym<coeff_t>(cpl, op, basis_in, fill);
    }
  } else {
    using bit_t = typename basis_t::bit_t;
    coeff_t mu = cpl.scalar().as<coeff_t>();
    int64_t s = op[0];
    auto apply = [&](bit_t ups, bit_t dns) {
      return bits::gbit(dns, s) ? mu : 0.;
    };
    generic_term_diag<symmetric, coeff_t>(basis_in, basis_out, apply, fill);
  }
}

} // namespace xdiag::basis::electron
