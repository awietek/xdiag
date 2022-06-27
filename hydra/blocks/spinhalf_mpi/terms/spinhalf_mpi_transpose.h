#pragma once

#include <vector>

#include <hydra/common.h>
#include <hydra/mpi/communicator.h>

namespace hydra::terms {

template <class bit_t>
mpi::Communicator spinhalf_mpi_transpose_com(SpinhalfMPI<bit_t> const &block,
                                             bool reverse = false) {
  using combinatorics::Combinations;

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_up = block.n_up();
  int n_postfix_bits = reverse ? block.n_prefix_bits_ : block.n_postfix_bits_;

  auto prefixes = reverse ? block.postfixes_ : block.prefixes_;
  auto postfixes = reverse ? block.prefixes_ : block.postfixes_;

  // Determine how many states I send
  std::vector<idx_t> n_states_i_send(mpi_size, 0);
  for (auto prefix : prefixes) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    for (auto postfix : Combinations<bit_t>(n_postfix_bits, n_up_postfix)) {
      int target_proc = block.process(postfix);
      ++n_states_i_send[target_proc];
    }
  }

  return mpi::Communicator(n_states_i_send);
}

template <class bit_t, class coeff_t>
void spinhalf_mpi_transpose(SpinhalfMPI<bit_t> const &block,
                            std::vector<coeff_t> const &vec,
                            std::vector<coeff_t> &send_buffer,
                            std::vector<coeff_t> &recv_buffer,
                            bool reverse = false) {
  using combinatorics::Combinations;
  using combinatorics::Subsets;

  int mpi_size, mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_up = block.n_up();
  int n_prefix_bits = reverse ? block.n_postfix_bits_ : block.n_prefix_bits_;
  int n_postfix_bits = reverse ? block.n_prefix_bits_ : block.n_postfix_bits_;

  auto prefixes = reverse ? block.postfixes_ : block.prefixes_;
  auto postfixes = reverse ? block.prefixes_ : block.postfixes_;

  auto &prefix_lintables =
      reverse ? block.postfix_lintables_ : block.prefix_lintables_;
  auto &postfix_limits = reverse ? block.prefix_limits_ : block.postfix_limits_;

  // Set up buffers for communicating
  mpi::Communicator com = spinhalf_mpi_transpose_com(block, reverse);
  idx_t buffer_size = std::max(com.recv_buffer_size(), com.send_buffer_size());
  send_buffer.resize(buffer_size);
  recv_buffer.resize(buffer_size);

  // Fill send buffer
  idx_t idx = 0;
  for (auto prefix : prefixes) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    for (auto postfix : Combinations<bit_t>(n_postfix_bits, n_up_postfix)) {
      int target_rank = block.process(postfix);
      com.add_to_send_buffer(target_rank, vec[idx], send_buffer);
      ++idx;
    }
  }

  // Communicate
  com.all_to_all(send_buffer, recv_buffer);

  // Sort reveived coefficients to postfix ordering (this is gnarly!!!)
  auto recv_offsets = com.n_values_i_recv_offsets();
  std::vector<idx_t> offsets(mpi_size, 0);
  for (auto prefix : Subsets<bit_t>(n_prefix_bits)) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits))
      continue;

    int origin_rank = block.process(prefix);
    idx_t origin_offset = recv_offsets[origin_rank];

    idx_t prefix_idx = prefix_lintables[n_up_prefix].index(prefix);
    for (auto postfix : postfixes) {
      if (bitops::popcnt(postfix) != n_up_postfix)
        continue;

      idx_t idx_received = origin_offset + offsets[origin_rank];
      idx_t idx_sorted = postfix_limits.at(postfix).first + prefix_idx;

      send_buffer[idx_sorted] = recv_buffer[idx_received];
      ++offsets[origin_rank];
    }
  }
  std::fill(recv_buffer.begin(), recv_buffer.end(), 0);
}

} // namespace hydra::terms
