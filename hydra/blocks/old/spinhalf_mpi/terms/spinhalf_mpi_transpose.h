#pragma once

#include <vector>

#include <hydra/common.h>
#include <hydra/mpi/buffer.h>
#include <hydra/mpi/communicator.h>

namespace hydra::terms {

template <class bit_t>
mpi::Communicator spinhalf_mpi_transpose_com(
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing,
    bool reverse = false) {

  // Determine how many states I send
  std::vector<int64_t> n_states_i_send(indexing.mpi_size(), 0);

  if (reverse) {
    for (auto postfix : indexing.postfixes()) {
      for (auto prefix : indexing.prefixes(postfix)) {
        int target_proc = indexing.process(prefix);
        ++n_states_i_send[target_proc];
      }
    }
  } else {

    for (auto prefix : indexing.prefixes()) {
      for (auto postfix : indexing.postfixes(prefix)) {
        int target_proc = indexing.process(postfix);
        ++n_states_i_send[target_proc];
      }
    }
  }

  return mpi::Communicator(n_states_i_send);
}

template <class bit_t, class coeff_t>
void spinhalf_mpi_transpose(
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing, const coeff_t *vec,
    bool reverse = false) {

  // Set up buffers for communicating
  mpi::Communicator com = spinhalf_mpi_transpose_com(indexing, reverse);
  int64_t buffer_size = std::max(com.recv_buffer_size(), com.send_buffer_size());
  mpi::buffer.reserve<coeff_t>(buffer_size);
  auto send_buffer = mpi::buffer.send<coeff_t>();
  auto recv_buffer = mpi::buffer.recv<coeff_t>();

  auto prefixes = reverse ? indexing.postfixes() : indexing.prefixes();

  // Fill send buffer
  int64_t idx = 0;
  for (auto prefix : prefixes) {
    auto postfixes =
        reverse ? indexing.prefixes(prefix) : indexing.postfixes(prefix);
    for (auto postfix : postfixes) {
      int target_rank = indexing.process(postfix);
      // lila::Log("bufsize {}", buffer_size);
      com.add_to_send_buffer(target_rank, vec[idx], send_buffer);
      ++idx;
    }
  }

  // Communicate
  com.all_to_all(send_buffer, recv_buffer);

  // Sort reveived coefficients to postfix ordering (this is gnarly!!!)
  auto recv_offsets = com.n_values_i_recv_offsets();
  std::vector<int64_t> offsets(indexing.mpi_size(), 0);

  int n_prefix_bits =
      reverse ? indexing.n_postfix_bits() : indexing.n_prefix_bits();
  int n_postfix_bits =
      reverse ? indexing.n_prefix_bits() : indexing.n_postfix_bits();

  for (auto prefix : combinatorics::Subsets<bit_t>(n_prefix_bits)) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = indexing.n_up() - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits))
      continue;

    int origin_rank = indexing.process(prefix);
    int64_t origin_offset = recv_offsets[origin_rank];
    int64_t prefix_idx = 0;
    if (reverse) {
      prefix_idx = indexing.postfix_indexing(((bit_t)1 << n_up_postfix) - 1)
                       .index(prefix);
    } else {
      prefix_idx = indexing.prefix_indexing(((bit_t)1 << n_up_postfix) - 1)
                       .index(prefix);
    }

    auto postfixes = reverse ? indexing.prefixes() : indexing.postfixes();
    for (auto postfix : postfixes) {
      if (bitops::popcnt(postfix) != n_up_postfix)
        continue;

      int64_t idx_received = origin_offset + offsets[origin_rank];
      int64_t postfix_begin = 0;
      if (reverse) {
        postfix_begin = indexing.prefix_begin(postfix);
      } else {
        postfix_begin = indexing.postfix_begin(postfix);
      }
      int64_t idx_sorted = postfix_begin + prefix_idx;

      send_buffer[idx_sorted] = recv_buffer[idx_received];
      ++offsets[origin_rank];
    }
  }
  mpi::buffer.clean_recv();
}

template <class bit_t, class coeff_t>
void spinhalf_mpi_transpose(
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing,
    std::vector<coeff_t> const &vec, bool reverse = false) {
  spinhalf_mpi_transpose(indexing, vec.data(), reverse);
}

} // namespace hydra::terms
