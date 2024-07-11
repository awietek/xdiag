#pragma once
#ifdef XDIAG_USE_MPI

#include <mpi.h>
#include <unordered_map>
#include <vector>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename bit_t, class process_f>
void fill_tables_sz(
    int64_t n_sites, int64_t nup, int64_t n_prefix_bits, process_f rank,
    std::vector<bit_t> &prefixes,
    std::unordered_map<bit_t, int64_t> &prefix_begin,
    std::vector<combinatorics::LinTable<bit_t>> &postfix_lintables,
    std::vector<std::vector<bit_t>> &postfix_states, ) {

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (nup < 0) { // Sz not conserved
    for (bit_t prefix : Subsets<bit_t>(n_prefix_bits)) {

      // only keep prefixes of this mpi rank
      if (rank(prefix) != mpi_rank) {
        continue;
      }
    }
  }

  // Determine the valid prefixes that belong to my process
  for (bit_t prefix : Subsets<bit_t>(n_prefix_bits)) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = n_up_ - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits_))
      continue;

    // only keep "random" prefixes
    if (process(prefix) != mpi_rank_)
      continue;

    prefix_begin_[prefix] = size_;
    size_ += combinatorics::binomial(n_postfix_bits_, n_up_postfix);
    prefixes_.push_back(prefix);
  }
}
} // namespace xdiag::basis::spinhalf_distributed

#endif
