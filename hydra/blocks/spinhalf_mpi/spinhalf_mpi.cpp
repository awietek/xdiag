#include "spinhalf_mpi.h"

#include <mpi.h>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/hashes.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/mpi/allreduce.h>

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/bitops.h>

namespace hydra {

template <class bit_t>
SpinhalfMPI<bit_t>::SpinhalfMPI(int n_sites, int n_up)
    : n_sites_(n_sites), n_prefix_bits_(n_sites / 2),
      n_postfix_bits_(n_sites - n_prefix_bits_), sz_conserved_(true),
      n_up_(n_up), n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_) {
  using combinatorics::Combinations;
  using combinatorics::Subsets;
  using namespace indexing;

  utils::check_nup_spinhalf(n_sites, n_up, "SpinhalfMPI");

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

  // Determine the valid prefixes that belong to my process
  size_ = 0;
  for (auto prefix : Subsets<bit_t>(n_prefix_bits_)) {
    int n_up_prefix = bitops::popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_postfix_bits_))
      continue;
    if (process(prefix) != mpi_rank_)
      continue;

    idx_t prefix_size = combinatorics::binomial(n_postfix_bits_, n_up_postfix);
    prefix_limits_[prefix] = {size_, size_ + prefix_size};
    size_ += prefix_size;
    prefixes_.push_back(prefix);
  }

  // Create the lintables for postfix lookup
  postfix_states_.resize(n_postfix_bits_ + 1);
  for (int n_up_postfix = 0; n_up_postfix <= n_postfix_bits_; ++n_up_postfix) {
    postfix_lintables_.push_back(
        LinTable<bit_t>(n_postfix_bits_, n_up_postfix));

    std::vector<bit_t> postfixes_nup(
        combinatorics::binomial(n_postfix_bits_, n_up_postfix), 0);
    idx_t pf_idx = 0;
    for (auto pf : Combinations<bit_t>(n_postfix_bits_, n_up_postfix)) {
      postfixes_nup[pf_idx++] = pf;
    }
    postfix_states_[n_up_postfix] = postfixes_nup;
  }

  // Determine the valid postfixes that belong to my process
  idx_t size = 0;
  for (auto postfix : Subsets<bit_t>(n_postfix_bits_)) {
    int n_up_postfix = bitops::popcnt(postfix);
    int n_up_prefix = n_up - n_up_postfix;
    if ((n_up_prefix < 0) || (n_up_prefix > n_prefix_bits_))
      continue;
    if (process(postfix) != mpi_rank_)
      continue;

    idx_t postfix_size = combinatorics::binomial(n_prefix_bits_, n_up_prefix);
    postfix_limits_[postfix] = {size, size + postfix_size};
    size += postfix_size;
    postfixes_.push_back(postfix);
  }

  // Create the lintables for postfix lookup
  prefix_states_.resize(n_prefix_bits_ + 1);
  for (int n_up_prefix = 0; n_up_prefix <= n_prefix_bits_; ++n_up_prefix) {
    prefix_lintables_.push_back(LinTable<bit_t>(n_prefix_bits_, n_up_prefix));

    std::vector<bit_t> prefixes_nup(
        combinatorics::binomial(n_prefix_bits_, n_up_prefix), 0);
    idx_t pf_idx = 0;
    for (auto pf : Combinations<bit_t>(n_prefix_bits_, n_up_prefix)) {
      prefixes_nup[pf_idx++] = pf;
    }
    prefix_states_[n_up_prefix] = prefixes_nup;
  }

  dim_ = 0;
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
