#pragma once

#include <algorithm>
#include <tuple>

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/models_mpi.h>
#include <hydra/models/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h>
#include <hydra/models/spinhalf_mpi/terms/spinhalf_mpi_exchange_mixed_com.h>
#include <hydra/models/spinhalf_mpi/terms/spinhalf_mpi_transpose.h>

#include <hydra/mpi/communicator.h>
#include <hydra/mpi/timing_mpi.h>

namespace hydra::terms::spinhalf_mpi {

template <class bit_t, class coeff_t>
void do_exchange_mpi(BondList const &bonds, Couplings const &couplings,
                     SpinhalfMPI<bit_t> const &block,
                     lila::Vector<coeff_t> const &vec_in,
                     lila::Vector<coeff_t> &vec_out) {
  using namespace indexing;
  using combinatorics::Combinations;
  using combinatorics::Subsets;

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  auto exchange_bonds = bonds.bonds_of_type("HEISENBERG") +
                        bonds.bonds_of_type("EXCHANGE") +
                        bonds.bonds_of_type("HB");

  int n_up = block.n_up();
  int n_prefix_bits = block.n_prefix_bits_;
  int n_postfix_bits = block.n_postfix_bits_;
  auto [prefix_bonds, postfix_bonds, mixed_bonds] =
      get_prefix_postfix_mixed_bonds(exchange_bonds, n_postfix_bits);

  /////////
  // Apply postfix bonds, no communication required
  auto tpost = rightnow_mpi();
  for (auto bond : postfix_bonds) {
    if (utils::coupling_is_zero(bond, couplings))
      continue;
    auto [s1, s2, Jhalf] = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
    assert(s1 < n_postfix_bits);
    assert(s2 < n_postfix_bits);
    bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    // loop through all postfixes
    idx_t idx = 0;
    idx_t idx_prefix = 0;
    for (auto prefix : block.prefixes_) {
      int n_up_prefix = utils::popcnt(prefix);
      int n_up_postfix = n_up - n_up_prefix;

      LinTable<bit_t> const &lintable = block.postfix_lintables_[n_up_postfix];
      auto const &postfixes = block.postfix_states_[n_up_postfix];
      for (auto postfix : postfixes) {

        if (utils::popcnt(postfix & mask) & 1) {
          bit_t new_postfix = postfix ^ mask;
          idx_t new_idx = idx_prefix + lintable.index(new_postfix);
          vec_out(new_idx) += Jhalf * vec_in(idx);
        }
        ++idx;
      }
      idx_prefix += combinatorics::binomial(n_postfix_bits, n_up_postfix);
      assert(idx_prefix == idx);
    }
  }
  timing_mpi(tpost, rightnow_mpi(), " (exchange) postfix", 2);

  /////////
  // Apply the prefix bonds after doing a transpose. then transpose back
  std::vector<coeff_t> send_buffer(vec_in.size(), 0);
  std::vector<coeff_t> recv_buffer(vec_in.size(), 0);
  auto tpre = rightnow_mpi();
  if (prefix_bonds.size() > 0) {

    // Transpose prefix/postfix order
    auto ttrans1 = rightnow_mpi();
    spinhalf_mpi_transpose(block, vec_in.vector(), send_buffer, recv_buffer,
                           false);
    timing_mpi(ttrans1, rightnow_mpi(), " (exchange) transpose 1", 2);

    // Loop over prefix_bonds
    for (auto bond : prefix_bonds) {
      if (utils::coupling_is_zero(bond, couplings))
        continue;
      auto [s1, s2, Jhalf] = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
      assert(s1 >= n_postfix_bits);
      assert(s2 >= n_postfix_bits);
      bit_t mask = ((bit_t)1 << (s1 - n_postfix_bits)) |
                   ((bit_t)1 << (s2 - n_postfix_bits));

      // loop through all prefixes
      idx_t idx = 0;
      idx_t idx_postfix = 0;
      for (auto postfix : block.postfixes_) {
        int n_up_postfix = utils::popcnt(postfix);
        int n_up_prefix = n_up - n_up_postfix;

        LinTable<bit_t> const &lintable = block.prefix_lintables_[n_up_prefix];
        auto const &prefixes = block.prefix_states_[n_up_prefix];
        for (auto prefix : prefixes) {

          if (utils::popcnt(prefix & mask) & 1) {
            bit_t new_prefix = prefix ^ mask;
            idx_t new_idx = idx_postfix + lintable.index(new_prefix);
            recv_buffer[new_idx] += Jhalf * send_buffer[idx];
          }
          ++idx;
        }
        idx_postfix += combinatorics::binomial(n_prefix_bits, n_up_prefix);
        assert(idx_postfix == idx);
      }
    }

    // Transpose back
    auto ttrans2 = rightnow_mpi();
    spinhalf_mpi_transpose(block, recv_buffer, send_buffer, recv_buffer, true);
    timing_mpi(ttrans2, rightnow_mpi(), " (exchange) transpose 2", 2);

    for (idx_t idx = 0; idx < vec_out.size(); ++idx) {
      vec_out(idx) += send_buffer[idx];
    }
  }

  timing_mpi(tpre, rightnow_mpi(), " (exchange) prefix", 2);

  /////////
  // Apply the mixed bonds
  auto tmix = rightnow_mpi();
  double tcom = 0;
  if (mixed_bonds.size() > 0) {

    // Figure out communication patterns and resize buffers accoringly
    auto [coms_mixed, max_send_size, max_recv_size] =
        spinhalf_mpi_exchange_mixed_com(block, mixed_bonds, couplings);

    if (max_send_size > send_buffer.size())
      send_buffer.resize(max_send_size);
    if (max_recv_size > recv_buffer.size())
      recv_buffer.resize(max_recv_size);

    int bond_idx = 0;
    for (auto bond : mixed_bonds) {

      if (utils::coupling_is_zero(bond, couplings))
        continue;
      auto [s1, s2, Jhalf] = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
      assert(s1 < n_postfix_bits);
      assert(s2 >= n_postfix_bits);
      bit_t prefix_mask = ((bit_t)1 << (s2 - n_postfix_bits));
      bit_t postfix_mask = ((bit_t)1 << s1);

      std::fill(send_buffer.begin(), send_buffer.end(), 0);
      std::fill(recv_buffer.begin(), recv_buffer.end(), 0);

      // Prepare sending states for this bond
      auto com = coms_mixed.at(bond_idx++);

      // Loop through all my states and fill them in send buffer
      idx_t idx = 0;
      for (auto prefix : block.prefixes_) {
        int n_up_prefix = utils::popcnt(prefix);
        int n_up_postfix = n_up - n_up_prefix;

        bit_t prefix_flipped = prefix ^ prefix_mask;
        int target_rank = block.process(prefix_flipped);

        auto const &postfixes = block.postfix_states_[n_up_postfix];

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
      double tcom1 = rightnow_mpi();
      com.all_to_all(send_buffer, recv_buffer);
      double tcom2 = rightnow_mpi();
      tcom += tcom2 - tcom1;

      // Fill received states into vec_out (gnarlyy!!!)
      // auto recv_offsets = com.n_values_i_recv_offsets();
      std::vector<idx_t> offsets(mpi_size, 0);
      for (auto prefix : Subsets<bit_t>(n_prefix_bits)) {

        // Only consider prefix if both itself and flipped version are valid
        int n_up_prefix = utils::popcnt(prefix);
        int n_up_postfix = n_up - n_up_prefix;
        if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits))
          continue;

        bit_t prefix_flipped = prefix ^ prefix_mask;
        int n_up_prefix_flipped = utils::popcnt(prefix_flipped);
        int n_up_postfix_flipped = n_up - n_up_prefix_flipped;
        if ((n_up_postfix_flipped < 0) ||
            (n_up_postfix_flipped > n_postfix_bits))
          continue;

        // Only consider prefix if it got sent to this mpi_rank
        int target_rank = block.process(prefix_flipped);
        if (target_rank != mpi_rank)
          continue;

        int origin_rank = block.process(prefix);
        idx_t origin_offset = recv_offsets[origin_rank];
        idx_t prefix_offset = block.prefix_limits_.at(prefix_flipped).first;

        auto const &postfixes = block.postfix_states_[n_up_postfix];

        // prefix up, postfix must be dn
        if (prefix & prefix_mask) {
          auto &postfix_lintable = block.postfix_lintables_[n_up_postfix + 1];

          for (auto postfix : postfixes) {
            if (!(postfix & postfix_mask)) {
              bit_t postfix_flipped = postfix ^ postfix_mask;
              idx_t idx_received = origin_offset + offsets[origin_rank];
              idx_t idx_target =
                  prefix_offset + postfix_lintable.index(postfix_flipped);
              vec_out(idx_target) += Jhalf * recv_buffer[idx_received];
              ++offsets[origin_rank];
            }
          }
        }
        // prefix dn, postfix must be up
        else {
          auto &postfix_lintable = block.postfix_lintables_[n_up_postfix - 1];
          for (auto postfix : postfixes) {
            if (postfix & postfix_mask) {
              bit_t postfix_flipped = postfix ^ postfix_mask;
              idx_t idx_received = origin_offset + offsets[origin_rank];
              idx_t idx_target =
                  prefix_offset + postfix_lintable.index(postfix_flipped);
              vec_out(idx_target) += Jhalf * recv_buffer[idx_received];
              ++offsets[origin_rank];
            }
          }
        }
      }
    } // for (auto bond : mixed_bonds)

    timing_mpi(tmix, rightnow_mpi(), " (exchange) mixed total", 2);
    LogMPI.out(2, " (exchange) mixed total (communication): {:.6f} secs", tcom);
  }
}

} // namespace hydra::terms::spinhalf_mpi
