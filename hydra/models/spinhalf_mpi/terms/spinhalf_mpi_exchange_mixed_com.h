#pragma once

#include <vector>

#include <hydra/common.h>
#include <hydra/mpi/communicator.h>
#include <hydra/models/models.h>
#include <hydra/models/utils/model_utils.h>
#include <hydra/models/utils/model_utils_mpi.h>


namespace hydra::spinhalfterms {

// Create a communicator for every mixed bond, and the maximal communicator size
template <class bit_t>
std::tuple<std::vector<mpi::Communicator>, idx_t, idx_t>
spinhalf_mpi_exchange_mixed_com(SpinhalfMPI<bit_t> const &block,
                                BondList const &mixed_bonds,
                                Couplings const &couplings) {

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_up = block.n_up();
  int n_postfix_bits = block.n_postfix_bits_;

  std::vector<mpi::Communicator> communicators;
  idx_t max_send_size = 0;
  idx_t max_recv_size = 0;
  for (auto bond : mixed_bonds) {

    std::vector<idx_t> n_states_i_send(mpi_size, 0);

    if (utils::coupling_is_zero(bond, couplings))
      continue;
    auto [s1, s2, Jhalf] = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
    assert(s1 < n_postfix_bits);
    assert(s2 >= n_postfix_bits);
    bit_t prefix_mask = ((bit_t)1 << (s2 - n_postfix_bits));
    bit_t postfix_mask = ((bit_t)1 << s1);

    for (auto prefix : block.prefixes_) {
      int n_up_prefix = utils::popcnt(prefix);
      int n_up_postfix = n_up - n_up_prefix;

      bit_t prefix_flipped = prefix ^ prefix_mask;
      int target_rank = block.process(prefix_flipped);

      // prefix up, postfix must be dn
      if (prefix & prefix_mask) {
        for (auto postfix : Combinations<bit_t>(n_postfix_bits, n_up_postfix)) {
          n_states_i_send[target_rank] += !(bool)(postfix & postfix_mask);
        }
      }

      // prefix dn, postfix must be dn
      if (!(prefix & prefix_mask)) {
        for (auto postfix : Combinations<bit_t>(n_postfix_bits, n_up_postfix)) {
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

} // namespace hydra::spinhalfterms
