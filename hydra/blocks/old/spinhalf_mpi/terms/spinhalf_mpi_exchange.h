#pragma once

#include <algorithm>
#include <tuple>

#include <hydra/indexing/spinhalf_mpi/spinhalf_mpi_indexing_sz.h>

#include <hydra/bitops/bitops.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

#include <hydra/blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h>
#include <hydra/blocks/spinhalf_mpi/terms/spinhalf_mpi_exchange_mixed_com.h>
#include <hydra/blocks/spinhalf_mpi/terms/spinhalf_mpi_transpose.h>

#include <hydra/mpi/communicator.h>
#include <hydra/mpi/timing_mpi.h>

namespace hydra::terms {

template <class bit_t, class coeff_t>
void spinhalf_mpi_exchange(
    BondList const &bonds, Couplings const &couplings,
    indexing::SpinhalfMPIIndexingSz<bit_t> const &indexing,
    lila::Vector<coeff_t> const &vec_in, lila::Vector<coeff_t> &vec_out) {
  using namespace indexing;
  using combinatorics::Combinations;
  using combinatorics::Subsets;

  assert(indexing.size() == vec_in.size());
  assert(indexing.size() == vec_out.size());

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  auto exchange_bonds = utils::clean_bondlist(
      bonds, couplings, {"HEISENBERG", "HB", "EXCHANGE"}, 2);
  int n_up = indexing.n_up();
  int n_prefix_bits = indexing.n_prefix_bits();
  int n_postfix_bits = indexing.n_postfix_bits();

  auto [prefix_bonds, postfix_bonds, mixed_bonds] =
      get_prefix_postfix_mixed_bonds(exchange_bonds, n_postfix_bits);

  /////////
  // Apply postfix bonds, no communication required
  auto tpost = rightnow_mpi();
  for (auto bond : postfix_bonds) {

    auto [s1, s2, Jhalf] = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
    assert(s1 < n_postfix_bits);
    assert(s2 < n_postfix_bits);
    bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

    // loop through all postfixes
    int64_t idx = 0;
    for (auto prefix : indexing.prefixes()) {
      auto const &lintable = indexing.postfix_indexing(prefix);
      auto const &postfixes = indexing.postfixes(prefix);
      int64_t prefix_begin = indexing.prefix_begin(prefix);
      for (auto postfix : postfixes) {

        if (bitops::popcnt(postfix & mask) & 1) {
          bit_t new_postfix = postfix ^ mask;
          int64_t new_idx = prefix_begin + lintable.index(new_postfix);
          vec_out(new_idx) += Jhalf * vec_in(idx);
        }
        ++idx;
      }
    }
  }
  timing_mpi(tpost, rightnow_mpi(), " (exchange) postfix", 2);

  /////////
  // Apply the prefix bonds after doing a transpose. then transpose back

  // prepare send/recv buffers
  int64_t buffer_size = indexing.size_max();
  mpi::buffer.reserve<coeff_t>(buffer_size);
  auto send_buffer = mpi::buffer.send<coeff_t>();
  auto recv_buffer = mpi::buffer.recv<coeff_t>();

  auto tpre = rightnow_mpi();
  if (prefix_bonds.size() > 0) {

    // Transpose prefix/postfix order
    auto ttrans1 = rightnow_mpi();
    spinhalf_mpi_transpose(indexing, vec_in.vector(), false);
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
      int64_t idx = 0;
      for (auto postfix : indexing.postfixes()) {
        auto const &lintable = indexing.prefix_indexing(postfix);
        auto const &prefixes = indexing.prefixes(postfix);
        int64_t postfix_begin = indexing.postfix_begin(postfix);
        for (auto prefix : prefixes) {

          if (bitops::popcnt(prefix & mask) & 1) {
            bit_t new_prefix = prefix ^ mask;
            int64_t new_idx = postfix_begin + lintable.index(new_prefix);
            recv_buffer[new_idx] += Jhalf * send_buffer[idx];
          }

          ++idx;
        }
      }
    }

    // Transpose back
    auto ttrans2 = rightnow_mpi();
    spinhalf_mpi_transpose(indexing, recv_buffer, true);
    timing_mpi(ttrans2, rightnow_mpi(), " (exchange) transpose 2", 2);

    for (int64_t idx = 0; idx < vec_out.size(); ++idx) {
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
        spinhalf_mpi_exchange_mixed_com(indexing, mixed_bonds, couplings);

    // prepare send/recv buffers
    mpi::buffer.reserve<coeff_t>(max_send_size, max_recv_size);
    auto send_buffer = mpi::buffer.send<coeff_t>();
    auto recv_buffer = mpi::buffer.recv<coeff_t>();

    int bond_idx = 0;
    for (auto bond : mixed_bonds) {

      if (utils::coupling_is_zero(bond, couplings))
        continue;
      auto [s1, s2, Jhalf] = get_exchange_s1_s2_Jhalf_mpi(bond, couplings);
      assert(s1 < n_postfix_bits);
      assert(s2 >= n_postfix_bits);
      bit_t prefix_mask = ((bit_t)1 << (s2 - n_postfix_bits));
      bit_t postfix_mask = ((bit_t)1 << s1);

      mpi::buffer.clean_send();
      mpi::buffer.clean_recv();
	
      // Prepare sending states for this bond
      auto com = coms_mixed.at(bond_idx++);

      // Loop through all my states and fill them in send buffer
      int64_t idx = 0;
      for (auto prefix : indexing.prefixes()) {
        bit_t prefix_flipped = prefix ^ prefix_mask;
        int target_rank = indexing.process(prefix_flipped);

        auto const &postfixes = indexing.postfixes(prefix);

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
      std::vector<int64_t> offsets(mpi_size, 0);
      for (auto prefix : Subsets<bit_t>(n_prefix_bits)) {

        // Only consider prefix if both itself and flipped version are valid
        int n_up_prefix = bitops::popcnt(prefix);
        int n_up_postfix = n_up - n_up_prefix;
        if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits))
          continue;

        bit_t prefix_flipped = prefix ^ prefix_mask;
        int n_up_prefix_flipped = bitops::popcnt(prefix_flipped);
        int n_up_postfix_flipped = n_up - n_up_prefix_flipped;
        if ((n_up_postfix_flipped < 0) ||
            (n_up_postfix_flipped > n_postfix_bits))
          continue;

        // Only consider prefix if it got sent to this mpi_rank
        int target_rank = indexing.process(prefix_flipped);
        if (target_rank != mpi_rank)
          continue;

        int origin_rank = indexing.process(prefix);
        int64_t origin_offset = recv_offsets[origin_rank];

        auto const &postfixes = indexing.postfixes(prefix);
        auto const &postfix_flipped_lintable =
            indexing.postfix_indexing(prefix_flipped);
        int64_t prefix_flipped_offset = indexing.prefix_begin(prefix_flipped);

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
      }
    } // for (auto bond : mixed_bonds)

    timing_mpi(tmix, rightnow_mpi(), " (exchange) mixed total", 2);
    LogMPI.out(2, " (exchange) mixed total (communication): {:.6f} secs", tcom);
  }
}

} // namespace hydra::terms
