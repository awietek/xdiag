#pragma once

#include <algorithm>
#include <tuple>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_exchange_postfix(Op const &op, basis_t const &basis,
                            arma::Col<coeff_t> const &vec_in,
                            arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis::bit_t;

  assert(basis.size() == vec_in.size());
  assert(basis.size() == vec_out.size());
  assert(op.type() == "EXCHANGE");
  assert(op.size() == 2);
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  assert((s1 >= 0) && (s2 >= 0));

  if (s1 == s2) {
    XDIAG_THROW("EXCHANGE Op with both sites equal not implemented yet");
  }

  assert(op.coupling_is<coeff_t>());
  coeff_t J = op.coupling<coeff_t>();
  coeff_t Jhalf = J / 2.0;

  int64_t n_prefix_bits = basis.n_prefix_bits();
  int64_t n_postfix_bits = basis.n_postfix_bits();
  assert((s1 < n_postfix_bits) && (s2 < n_postfix_bits));

  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  int64_t idx = 0;
  for (bit_t prefix : basis.prefixes()) {
    auto const &lintable = basis.postfix_indexing(prefix);
    auto const &postfixes = basis.postfixes(prefix);
    int64_t prefix_begin = basis.prefix_begin(prefix);
    for (bit_t postfix : postfixes) {

      if (bitops::popcnt(postfix & mask) & 1) {
        bit_t new_postfix = postfix ^ mask;
        int64_t new_idx = prefix_begin + lintable.index(new_postfix);
        vec_out(new_idx) += Jhalf * vec_in(idx);
      }
      ++idx;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <class basis_t, typename coeff_t>
void apply_exchange_prefix(Op const &op, basis_t const &basis,
                           arma::Col<coeff_t> const &vec_in,
                           arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis::bit_t;

  assert(basis.size() == vec_in.size());
  assert(basis.size() == vec_out.size());
  assert(op.type() == "EXCHANGE");
  assert(op.size() == 2);
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  assert((s1 >= 0) && (s2 >= 0));

  if (s1 == s2) {
    XDIAG_THROW("EXCHANGE Op with both sites equal not implemented yet");
  }

  assert(op.coupling_is<coeff_t>());
  coeff_t J = op.coupling<coeff_t>();
  coeff_t Jhalf = J / 2.0;

  int64_t n_prefix_bits = basis.n_prefix_bits();
  int64_t n_postfix_bits = basis.n_postfix_bits();
  assert((s1 >= n_postfix_bits) && (s2 += n_postfix_bits));

  bit_t mask = (((bit_t)1 << s1) | ((bit_t)1 << s2)) >> n_postfix_bits;

  // loop through all postfixes
  int64_t idx = 0;
  for (bit_t postfix : basis.postfixes()) {
    auto const &lintable = basis.prefix_indexing(postfix);
    auto const &prefixes = basis.prefixes(postfix);
    int64_t postfix_begin = basis.postfix_begin(postfix);
    for (bit prefix : prefixes) {
      if (bitops::popcnt(prefix & mask) & 1) {
        bit_t new_prefix = prefix ^ mask;
        int64_t new_idx = postfix_begin + lintable.index(new_prefix);
        recv_buffer[new_idx] += Jhalf * send_buffer[idx];
      }
      ++idx;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <class basis_t, typename coeff_t>
void apply_exchange_mixed(Op const &op, basis_t const &basis,
                          arma::Col<coeff_t> const &vec_in,
                          arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis::bit_t;

  int32_t mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  assert(basis.size() == vec_in.size());
  assert(basis.size() == vec_out.size());
  assert(op.type() == "EXCHANGE");
  assert(op.size() == 2);
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  assert((s1 >= 0) && (s2 >= 0));

  if (s1 == s2) {
    XDIAG_THROW("EXCHANGE Op with both sites equal not implemented yet");
  }
  int64_t ss1 = std::min(s1, s2);
  int64_t ss2 = std::max(s1, s2);
  bit_t prefix_mask = ((bit_t)1 << (ss2 - n_postfix_bits));
  bit_t postfix_mask = ((bit_t)1 << s11);

  assert(op.coupling_is<coeff_t>());
  coeff_t J = op.coupling<coeff_t>();
  coeff_t Jhalf = J / 2.0;

  int64_t n_postfix_bits = basis.n_postfix_bits();

  mpi::Communicator comm;
  // Check whether communication pattern has already been determined
  if (basis.comm_pattern().contains(op)) {
    comm = basis.comm_pattern[op];
  }
  // if not, compute it anew
  else {
    std::vector<int64_t> n_states_i_send(mpi_size, 0);

    for (bit_t prefix : basis.prefixes()) {

      bit_t prefix_flipped = prefix ^ prefix_mask;
      int32_t target_rank = basis.rank(prefix_flipped);

      // prefix up, postfix must be dn
      if (prefix & prefix_mask) {
        for (bit_t postfix : basis.postfixes(prefix)) {
          n_states_i_send[target_rank] += !(bool)(postfix & postfix_mask);
        }
      }

      // prefix dn, postfix must be dn
      if (!(prefix & prefix_mask)) {
        for (auto postfix : basis.postfixes(prefix)) {
          n_states_i_send[target_rank] += (bool)(postfix & postfix_mask);
        }
      }
    }

    auto comm = mpi::Communicator(n_states_i_send);
    basis.comm_pattern[op] = comm;
  }

  // prepare send/recv buffers
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
    auto const &postfixes = basis.postfixes(prefix);

    // prefix up, postfix must be dn
    if (prefix & prefix_mask) {
      for (auto postfix : postfixes) {
        if (!(postfix & postfix_mask)) {
          com.add_to_send_buffer(target_rank, vec_in(idx), send_buffer);
        }
        ++idx;
      }
    }
    // prefix dn, postfix must be up
    else {
      for (auto postfix : postfixes) {
        if (postfix & postfix_mask) {
          com.add_to_send_buffer(target_rank, vec_in(idx), send_buffer);
        }
        ++idx;
      }
    }
  }
  auto send_offsets = com.n_values_i_send_offsets();
  auto recv_offsets = com.n_values_i_recv_offsets();

  // Communicate
  com.all_to_all(send_buffer, recv_buffer);

  // Fill received states into vec_out (gnarlyy!!!)
  // auto recv_offsets = com.n_values_i_recv_offsets();
  std::vector<int64_t> offsets(mpi_size, 0);
  for (bit_t prefix : Subsets<bit_t>(n_prefix_bits)) {

    // Only consider prefix if both itself and flipped version are valid
    int64_t n_up_prefix = bitops::popcnt(prefix);
    int64_t n_up_postfix = n_up - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits)) {
      continue;
    }

    bit_t prefix_flipped = prefix ^ prefix_mask;
    int64_t n_up_prefix_flipped = bitops::popcnt(prefix_flipped);
    int64_t n_up_postfix_flipped = n_up - n_up_prefix_flipped;
    if ((n_up_postfix_flipped < 0) || (n_up_postfix_flipped > n_postfix_bits))
      continue;

    // Only consider prefix if it got sent to this mpi_rank
    int target_rank = basis.rank(prefix_flipped);
    if (target_rank != mpi_rank)
      continue;

    int32_t origin_rank = basis.rank(prefix);
    int64_t origin_offset = recv_offsets[origin_rank];

    auto const &postfixes = basis.postfixes(prefix);
    auto const &postfix_flipped_lintable =
        basis.postfix_indexing(prefix_flipped);
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
          vec_out(int64_target) += Jhalf * recv_buffer[idx_received];
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
          vec_out(int64_target) += Jhalf * recv_buffer[idx_received];
          ++offsets[origin_rank];
        }
      }
    }
  } // for (bit_t prefix : Subsets<bit_t>(n_prefix_bits))
}

} // namespace xdiag::basis::spinhalf_distributed
