// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "transpose.hpp"

#include <xdiag/combinatorics/subsets.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class bit_t, typename coeff_t>
void transpose(BasisSz<bit_t> const &basis, coeff_t const *vec_in,
               bool reverse) {
  mpi::Communicator com = basis.transpose_communicator(reverse);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Adjust the global MPI buffer size if necessary
  int64_t buffer_size =
      std::max(com.recv_buffer_size(), com.send_buffer_size());
  mpi::buffer.reserve<coeff_t>(buffer_size);

  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  auto const& prefixes = reverse ? basis.postfixes() : basis.prefixes();

  // Fill send buffer
  int64_t idx = 0;
  for (auto prefix : prefixes) {
    auto const& postfixes = reverse ? basis.prefix_states(prefix) : basis.postfix_states(prefix);
    for (auto postfix : postfixes) {
      int target_rank = basis.rank(postfix);
      com.add_to_send_buffer(target_rank, vec_in[idx], send_buffer);
      ++idx;
    }
  }

  // Communicate
  com.all_to_all(send_buffer, recv_buffer);

  // Sort reveived coefficients to postfix ordering (this is gnarly!!!)
  std::vector<int64_t> offsets(mpi_size, 0);

  int n_prefix_bits = reverse ? basis.n_postfix_bits() : basis.n_prefix_bits();
  int n_postfix_bits = reverse ? basis.n_prefix_bits() : basis.n_postfix_bits();

  for (auto prefix : combinatorics::Subsets<bit_t>(n_prefix_bits)) {
    int nup_prefix = bits::popcnt(prefix);
    int nup_postfix = basis.nup() - nup_prefix;
    if ((nup_postfix < 0) || (nup_postfix > n_postfix_bits))
      continue;

    int origin_rank = basis.rank(prefix);
    int64_t origin_offset = com.n_values_i_recv_offset(origin_rank);
    int64_t prefix_idx = 0;
    if (reverse) {
      bit_t postfix = ((bit_t)1 << nup_postfix) - 1;
      prefix_idx = basis.postfix_lintable(postfix).index(prefix);
    } else {
      bit_t postfix = ((bit_t)1 << nup_postfix) - 1;
      prefix_idx = basis.prefix_lintable(postfix).index(prefix);
    }

    auto postfixes = reverse ? basis.prefixes() : basis.postfixes();
    for (bit_t postfix : postfixes) {
      if (bits::popcnt(postfix) != nup_postfix)
        continue;

      int64_t idx_received = origin_offset + offsets[origin_rank];
      int64_t postfix_begin = 0;
      if (reverse) {
        postfix_begin = basis.prefix_begin(postfix);
      } else {
        postfix_begin = basis.postfix_begin(postfix);
      }
      int64_t idx_sorted = postfix_begin + prefix_idx;

      send_buffer[idx_sorted] = recv_buffer[idx_received];
      ++offsets[origin_rank];
    }
  }
  mpi::buffer.clean_recv();
}

template void transpose(BasisSz<uint32_t> const &, double const *, bool);
template void transpose(BasisSz<uint64_t> const &, double const *, bool);
template void transpose(BasisSz<uint32_t> const &, complex const *, bool);
template void transpose(BasisSz<uint64_t> const &, complex const *, bool);

} // namespace xdiag::basis::spinhalf_distributed
