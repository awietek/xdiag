// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/op.hpp>
#include <xdiag/parallel/mpi/buffer.hpp>
namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_spsm_postfix(Coupling const &cpl, Op const &op,
                        basis_t const &basis_in,
                        arma::Col<coeff_t> const &vec_in,
                        basis_t const &basis_out,
                        arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t H = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  std::string type = op.type();

  bit_t mask = ((bit_t)1 << s);

  int64_t n_postfix_bits = basis_in.n_postfix_bits();
  int64_t idx = 0;
  for (auto prefix : basis_in.prefixes()) {
    auto const &postfixes = basis_in.postfix_states(prefix);
    auto const &lintable = basis_out.postfix_lintable(prefix);
    int64_t idx_prefix = basis_out.prefix_begin(prefix);

    if (idx_prefix != invalid_index) { // can happen since total Sz changes
      if (type == "S+") {
        for (auto postfix_in : postfixes) {
          if (!(postfix_in & mask)) {
            bit_t postfix_out = postfix_in | mask;
            int64_t idx_out = idx_prefix + lintable.index(postfix_out);
            vec_out(idx_out) += H * vec_in(idx);
          }
          ++idx;
        }
      } else if (type == "S-") {
        for (auto postfix_in : postfixes) {
          if (postfix_in & mask) {
            bit_t postfix_out = postfix_in ^ mask;
            int64_t idx_out = idx_prefix + lintable.index(postfix_out);
            vec_out(idx_out) += H * vec_in(idx);
          }
          ++idx;
        }
      }
    } else {
      idx += postfixes.size();
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <class basis_t, typename coeff_t>
void apply_spsm_prefix(Coupling const &cpl, Op const &op,
                       basis_t const &basis_in, basis_t const &basis_out) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t H = cpl.scalar().as<coeff_t>();
  std::string type = op.type();
  int64_t s = op[0];

  int64_t n_prefix_bits = basis_in.n_prefix_bits();
  int64_t n_postfix_bits = basis_in.n_postfix_bits();
  bit_t mask = ((bit_t)1 << (s - n_postfix_bits));

  int64_t buffer_size = std::max(basis_out.size_max(), basis_in.size_max());
  mpi::buffer.reserve<coeff_t>(buffer_size);
  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  // loop through all postfixes
  int64_t idx = 0;
  for (auto postfix : basis_in.postfixes()) {

    auto const &prefixes = basis_in.prefix_states(postfix);
    auto const &lintable = basis_out.prefix_lintable(postfix);
    int64_t idx_postfix = basis_out.postfix_begin(postfix);

    if (idx_postfix != invalid_index) { // can happen since total Sz changes
      if (type == "S+") {
        for (auto prefix_in : prefixes) {
          if (!(prefix_in & mask)) {
            bit_t prefix_out = prefix_in | mask;
            int64_t idx_out = idx_postfix + lintable.index(prefix_out);
            recv_buffer[idx_out] += H * send_buffer[idx];
          }
          ++idx;
        }

      } else if (type == "S-") {
        for (auto prefix_in : prefixes) {
          if (prefix_in & mask) {
            bit_t prefix_out = prefix_in ^ mask;
            int64_t idx_out = idx_postfix + lintable.index(prefix_out);
            recv_buffer[idx_out] += H * send_buffer[idx];
          }
          ++idx;
        }
      }
    } else {
      idx += prefixes.size();
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag::basis::spinhalf_distributed
