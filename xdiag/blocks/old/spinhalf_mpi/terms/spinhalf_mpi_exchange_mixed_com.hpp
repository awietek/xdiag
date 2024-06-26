#pragma once

#include <vector>
#include <tuple>

#include <xdiag/common.hpp>
#include <xdiag/mpi/communicator.hpp>

#include <xdiag/indexing/spinhalf_mpi/spinhalf_mpi_indexing_sz.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/utils/block_utils.hpp>
#include <xdiag/blocks/utils/block_utils_mpi.hpp>

namespace xdiag::terms {

// Create a communicator for every mixed bond, and the maximal communicator size
template <class bit_t>
std::tuple<std::vector<mpi::Communicator>, int64_t, int64_t>
spinhalf_mpi_exchange_mixed_com(
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing,
    BondList const &mixed_bonds, Couplings const &couplings) {
  using combinatorics::Combinations;

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_postfix_bits = indexing.n_postfix_bits();

  std::vector<mpi::Communicator> communicators;
  int64_t max_send_size = 0;
  int64_t max_recv_size = 0;
  for (auto bond : mixed_bonds) {

    std::vector<int64_t> n_states_i_send(mpi_size, 0);

    int s1, s2;
    std::tie(s1, s2, std::ignore) = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
    assert(s1 < n_postfix_bits);
    assert(s2 >= n_postfix_bits);
    bit_t prefix_mask = ((bit_t)1 << (s2 - n_postfix_bits));
    bit_t postfix_mask = ((bit_t)1 << s1);

    for (auto prefix : indexing.prefixes()) {

      bit_t prefix_flipped = prefix ^ prefix_mask;
      int target_rank = indexing.process(prefix_flipped);

      // prefix up, postfix must be dn
      if (prefix & prefix_mask) {
        for (auto postfix : indexing.postfixes(prefix)) {
          n_states_i_send[target_rank] += !(bool)(postfix & postfix_mask);
        }
      }

      // prefix dn, postfix must be dn
      if (!(prefix & prefix_mask)) {
        for (auto postfix : indexing.postfixes(prefix)) {
          n_states_i_send[target_rank] += (bool)(postfix & postfix_mask);
        }
      }
    }

    auto com = mpi::Communicator(n_states_i_send);
    communicators.push_back(com);
    if (com.send_buffer_size() > max_send_size)
      max_send_size = com.send_buffer_size();
    if (com.recv_buffer_size() > max_recv_size)
      max_recv_size = com.recv_buffer_size();
  }
  return {communicators, max_send_size, max_recv_size};
}

} // namespace xdiag::terms
