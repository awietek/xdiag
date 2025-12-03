// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "norm.hpp"

#include <cmath>

namespace xdiag::symmetries {

template <typename bit_t, typename coeff_t, class group_action_t>
double norm(bit_t state, group_action_t const &group_action,
            std::vector<coeff_t> const &characters) {
  coeff_t amplitude = 0.0;
  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      amplitude += characters(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template <typename bit_t, typename coeff_t, class group_action_t>
double norm_fermionic(bit_t state, group_action_t const &group_action,
                      std::vector<coeff_t> const &characters) {
  coeff_t amplitude = 0.0;
  int64_t nsites = group_action.nsites();
  auto work = fermi_work(nsites);
  auto const &group = group_action.permutation_group();

  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {

    auto const &perm = group[sym];

    bit_t tstate = group_action.apply(sym, state);
    if (tstate == state) {
      if (fermi_bool_of_permutation(state, perm, work)) {
        amplitude -= characters(sym);
      } else {
        amplitude += characters(sym);
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template <typename bit_t, typename coeff_t, class group_action_t>
double norm_electron(bit_t ups, bit_t dns, group_action_t const &group_action,
                     std::vector<coeff_t> const &characters) {
  coeff_t amplitude = 0.0;
  int64_t nsites = group_action.nsites();
  auto work = fermi_work(nsites);
  auto const &group = group_action.permutation_group();

  for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {

    auto const &perm = group[sym];
    bit_t tups = group_action.apply(sym, ups);

    if (tups == ups) {
      bool fermi_bool_ups = fermi_bool_of_permutation(ups, perm, work);

      bit_t tdns = group_action.apply(sym, dns);

      if (tdns == dns) {
        bool fermi_bool_dns = fermi_bool_of_permutation(dns, perm, work);

        if (fermi_bool_ups == fermi_bool_dns) {
          amplitude += characters(sym);
        } else {
          amplitude -= characters(sym);
        }
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

template <typename bit_t, typename coeff_t, class group_action_t>
double norm_electron_subset(bit_t ups, bit_t dns,
                            group_action_t const &group_action,
                            std::vector<coeff_t> const &characters,
                            gsl::span<int64_t const> syms) {
  coeff_t amplitude = 0.0;
  int64_t nsites = group_action.nsites();
  auto work = fermi_work(nsites);
  auto const &group = group_action.permutation_group();

  for (int64_t sym : syms) {
    auto const &perm = group[sym];
    bit_t tups = group_action.apply(sym, ups);

    if (tups == ups) {
      bool fermi_bool_ups = fermi_bool_of_permutation(ups, perm, work);

      bit_t tdns = group_action.apply(sym, dns);

      if (tdns == dns) {
        bool fermi_bool_dns = fermi_bool_of_permutation(dns, perm, work);
        if (fermi_bool_ups == fermi_bool_dns) {
          amplitude += characters(sym);
        } else {
          amplitude -= characters(sym);
        }
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

} // namespace xdiag::symmetries
