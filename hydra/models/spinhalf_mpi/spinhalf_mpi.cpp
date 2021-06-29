#include "spinhalf_mpi.h"

#include <mpi.h>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/hashes.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/mpi/allreduce.h>

#include <hydra/utils/bitops.h>

namespace hydra {

template <class bit_t>
SpinhalfMPI<bit_t>::SpinhalfMPI(int n_sites, int n_up)
    : n_sites_(n_sites), n_prefix_bits_(n_sites / 2),
      n_postfix_bits_(n_sites - n_prefix_bits_), sz_conserved_(true),
      n_up_(n_up), n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_) {

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

  // Determine the valid prefixes that belong to my process
  idx_t idx = 0;
  for (auto prefix : Subsets<bit_t>(n_prefix_bits_)) {
    int n_up_prefix = utils::popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    if (n_up_postfix < 0)
      continue;
    if (process(prefix) != mpi_rank_)
      continue;

    idx_t start_idx = idx;
    idx += combinatorics::binomial(n_postfix_bits_, n_up_postfix);
    idx_t end_idx = idx;
    prefixes_.push_back(prefix);
    prefix_limits_[prefix] = {start_idx, end_idx};
  }

  // Create the lintables for postfix lookup
  for (int n_up_postfix = 0; n_up_postfix <= n_postfix_bits_; ++n_up_postfix) {
    postfix_lintables_.push_back(
        LinTable<bit_t>(n_postfix_bits_, n_up_postfix));
  }

  size_ = idx;
  mpi::Allreduce(&size_, &dim_, 1, MPI_SUM, MPI_COMM_WORLD);
}

template <class bit_t>
bool SpinhalfMPI<bit_t>::operator==(SpinhalfMPI<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <class bit_t>
bool SpinhalfMPI<bit_t>::operator!=(SpinhalfMPI<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class SpinhalfMPI<uint16>;
template class SpinhalfMPI<uint32>;
template class SpinhalfMPI<uint64>;

} // namespace hydra
