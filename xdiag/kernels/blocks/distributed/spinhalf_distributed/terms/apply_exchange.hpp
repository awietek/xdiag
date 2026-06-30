// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <tuple>

#include <xdiag/armadillo.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/mpi/buffer.hpp>
#include <xdiag/mpi/comm_pattern.hpp>
#include <xdiag/mpi/communicator.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_exchange_postfix(Coeff const &cpl, Op const &op,
                            basis_t const &basis,
                            arma::Col<coeff_t> const &vec_in,
                            arma::Col<coeff_t> &vec_out) {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  coeff_t Jhalf = J / 2.0;
  bool is_asym = (op.type() == "ExchangeAsym");
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  int64_t n_prefix_bits = basis.n_prefix_bits();
  int64_t n_postfix_bits = basis.n_postfix_bits();
  assert((s1 < n_postfix_bits) && (s2 < n_postfix_bits));

  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  bit_t s1mask = (bit_t)1 << s1;

  // Per-branch coupling, hoisted out of the loop. ExchangeAsym =
  // 1/2(S+_{s1}S-_{s2} - S-_{s1}S+_{s2}): the S-_{s1} branch (input s1 up) gets
  // -Jhalf (matching term_exchange_asym); plain Exchange uses Jhalf for both.
  coeff_t j_s1dn = Jhalf;
  coeff_t j_s1up = is_asym ? -Jhalf : Jhalf;

  int64_t idx = 0;
  for (bit_t prefix : basis.prefixes()) {
    auto const &lintable = basis.postfix_lintable(prefix);
    auto const &postfixes = basis.postfix_states(prefix);
    int64_t prefix_begin = basis.prefix_begin(prefix);
    for (bit_t postfix : postfixes) {

      if (bits::popcount(postfix & mask) & 1) {
        bit_t new_postfix = postfix ^ mask;
        int64_t new_idx = prefix_begin + lintable.index(new_postfix);
        vec_out(new_idx) += ((postfix & s1mask) ? j_s1up : j_s1dn) * vec_in(idx);
      }
      ++idx;
    }
  }
}

// apply_exchange_prefix works with the send/recv buffer as input/output vectors
template <class basis_t, typename coeff_t>
void apply_exchange_prefix(Coeff const &cpl, Op const &op,
                           basis_t const &basis) {
  using bit_t = typename basis_t::bit_t;

  int64_t buffer_size = basis.size_max();
  mpi::buffer.reserve<coeff_t>(buffer_size);
  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  coeff_t J = cpl.scalar().as<coeff_t>();
  coeff_t Jhalf = J / 2.0;
  bool is_asym = (op.type() == "ExchangeAsym");
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  int64_t n_prefix_bits = basis.n_prefix_bits();
  int64_t n_postfix_bits = basis.n_postfix_bits();
  assert((s1 >= n_postfix_bits) && (s2 >= n_postfix_bits));

  bit_t mask = (((bit_t)1 << s1) | ((bit_t)1 << s2)) >> n_postfix_bits;
  bit_t s1mask = (bit_t)1 << (s1 - n_postfix_bits);

  // Per-branch coupling, hoisted out of the loop (see apply_exchange_postfix):
  // ExchangeAsym flips the sign of the input-s1-up branch, plain Exchange does
  // not.
  coeff_t j_s1dn = Jhalf;
  coeff_t j_s1up = is_asym ? -Jhalf : Jhalf;

  // loop through all postfixes
  int64_t idx = 0;
  for (bit_t postfix : basis.postfixes()) {
    auto const &lintable = basis.prefix_lintable(postfix);
    auto const &prefixes = basis.prefix_states(postfix);
    int64_t postfix_begin = basis.postfix_begin(postfix);
    for (bit_t prefix : prefixes) {
      if (bits::popcount(prefix & mask) & 1) {
        bit_t new_prefix = prefix ^ mask;
        int64_t new_idx = postfix_begin + lintable.index(new_prefix);
        recv_buffer[new_idx] +=
            ((prefix & s1mask) ? j_s1up : j_s1dn) * send_buffer[idx];
      }
      ++idx;
    }
  }
}

template <class basis_t, typename coeff_t>
void apply_exchange_mixed(Coeff const &cpl, Op const &op,
                          basis_t const &basis,
                          arma::Col<coeff_t> const &vec_in,
                          arma::Col<coeff_t> &vec_out) {
  using bit_t = typename basis_t::bit_t;
  int32_t mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  coeff_t J = cpl.scalar().as<coeff_t>();
  coeff_t Jhalf = J / 2.0;
  bool is_asym = (op.type() == "ExchangeAsym");
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  assert((s1 >= 0) && (s2 >= 0));

  int64_t nup = basis.nup();
  int64_t n_prefix_bits = basis.n_prefix_bits();
  int64_t n_postfix_bits = basis.n_postfix_bits();
  int64_t ss1 = std::min(s1, s2);
  int64_t ss2 = std::max(s1, s2);
  bit_t prefix_mask = ((bit_t)1 << (ss2 - n_postfix_bits));
  bit_t postfix_mask = ((bit_t)1 << ss1);

  // Per-branch coupling, hoisted out of the fill loop. ExchangeAsym =
  // 1/2(S+_{s1}S-_{s2} - S-_{s1}S+_{s2}) gives -Jhalf when the input s1 spin is
  // up (matching term_exchange_asym); plain Exchange uses Jhalf throughout.
  // s1 lives on the prefix iff it is the larger site (ss2); the two fill
  // branches below fix the s1 spin from the prefix/postfix occupation.
  coeff_t j_s1up = is_asym ? -Jhalf : Jhalf;
  bool s1_on_prefix = (s1 > s2);
  // "prefix up, postfix dn" branch: s1 is up iff s1 sits on the prefix.
  coeff_t coeff_prefix_up = s1_on_prefix ? j_s1up : Jhalf;
  // "prefix dn, postfix up" branch: s1 is up iff s1 sits on the postfix.
  coeff_t coeff_prefix_dn = s1_on_prefix ? Jhalf : j_s1up;

  mpi::Communicator comm;

  // Check whether communication pattern has already been determined
  if (basis.comm_pattern().contains(op)) {
    comm = basis.comm_pattern()[op];
  }
  // if not, compute it anew
  else {
    std::vector<int64_t> n_states_i_send(mpi_size, 0);

    for (bit_t prefix : basis.prefixes()) {

      bit_t prefix_flipped = prefix ^ prefix_mask;
      int32_t target_rank = basis.rank(prefix_flipped);

      auto const &postfixes = basis.postfix_states(prefix);

      // prefix up, postfix must be dn
      if (prefix & prefix_mask) {
        for (bit_t postfix : postfixes) {
          n_states_i_send[target_rank] += !(bool)(postfix & postfix_mask);
        }
      }

      // prefix dn, postfix must be dn
      if (!(prefix & prefix_mask)) {
        for (auto postfix : postfixes) {
          n_states_i_send[target_rank] += (bool)(postfix & postfix_mask);
        }
      }
    }
    comm = mpi::Communicator(n_states_i_send);
    basis.comm_pattern().append(op, comm);
  }

  // prepare send/recv buffers
  int64_t max_send_size = comm.send_buffer_size();
  int64_t max_recv_size = comm.recv_buffer_size();
  mpi::buffer.reserve<coeff_t>(max_send_size, max_recv_size);
  auto send_buffer = mpi::buffer.send<coeff_t>();
  auto recv_buffer = mpi::buffer.recv<coeff_t>();
  mpi::buffer.clean_send();
  mpi::buffer.clean_recv();

  // Loop through all my states and fill them in send buffer
  int64_t idx = 0;
  for (auto prefix : basis.prefixes()) {
    bit_t prefix_flipped = prefix ^ prefix_mask;
    int32_t target_rank = basis.rank(prefix_flipped);
    auto const &postfixes = basis.postfix_states(prefix);

    // prefix up, postfix must be dn
    if (prefix & prefix_mask) {
      for (auto postfix : postfixes) {
        if (!(postfix & postfix_mask)) {
          comm.add_to_send_buffer(target_rank, vec_in(idx), send_buffer);
        }
        ++idx;
      }
    }
    // prefix dn, postfix must be up
    else {
      for (auto postfix : postfixes) {
        if (postfix & postfix_mask) {
          comm.add_to_send_buffer(target_rank, vec_in(idx), send_buffer);
        }
        ++idx;
      }
    }
  }

  // Communicate
  comm.all_to_all(send_buffer, recv_buffer);

  // Fill received states into vec_out (gnarlyy!!!)
  // auto recv_offsets = comm.n_values_i_recv_offsets();
  std::vector<int64_t> offsets(mpi_size, 0);
  for (bit_t prefix : combinatorics::Subsets<bit_t>(n_prefix_bits)) {

    // Only consider prefix if both itself and flipped version are valid
    int64_t nup_prefix = bits::popcount(prefix);
    int64_t nup_postfix = nup - nup_prefix;
    if ((nup_postfix < 0) || (nup_postfix > n_postfix_bits)) {
      continue;
    }

    bit_t prefix_flipped = prefix ^ prefix_mask;
    int64_t nup_prefix_flipped = bits::popcount(prefix_flipped);
    int64_t nup_postfix_flipped = nup - nup_prefix_flipped;
    if ((nup_postfix_flipped < 0) || (nup_postfix_flipped > n_postfix_bits))
      continue;

    // Only consider prefix if it got sent to this mpi_rank
    int target_rank = basis.rank(prefix_flipped);
    if (target_rank != mpi_rank)
      continue;

    int32_t origin_rank = basis.rank(prefix);
    int64_t origin_offset = comm.n_values_i_recv_offset(origin_rank);

    auto const &postfixes = basis.postfix_states(prefix);
    auto const &postfix_flipped_lintable =
        basis.postfix_lintable(prefix_flipped);
    int64_t prefix_flipped_offset = basis.prefix_begin(prefix_flipped);

    // prefix up, postfix must be dn
    if (prefix & prefix_mask) {
      for (auto postfix : postfixes) {
        if (!(postfix & postfix_mask)) {
          bit_t postfix_flipped = postfix ^ postfix_mask;
          int64_t idx_received = origin_offset + offsets[origin_rank];
          int64_t int64_target =
              prefix_flipped_offset +
              postfix_flipped_lintable.index(postfix_flipped);
          vec_out(int64_target) += coeff_prefix_up * recv_buffer[idx_received];
          ++offsets[origin_rank];
        }
      }
    }
    // prefix dn, postfix must be up
    else {
      for (auto postfix : postfixes) {
        if (postfix & postfix_mask) {
          bit_t postfix_flipped = postfix ^ postfix_mask;
          int64_t idx_received = origin_offset + offsets[origin_rank];
          int64_t int64_target =
              prefix_flipped_offset +
              postfix_flipped_lintable.index(postfix_flipped);
          vec_out(int64_target) += coeff_prefix_dn * recv_buffer[idx_received];
          ++offsets[origin_rank];
        }
      }
    }
  } // for (bit_t prefix : Subsets<bit_t>(n_prefix_bits))
}

} // namespace xdiag::basis::spinhalf_distributed
