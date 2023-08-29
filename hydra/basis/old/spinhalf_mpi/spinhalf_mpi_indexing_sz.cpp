#include "spinhalf_mpi_indexing_sz.h"

#include <mpi.h>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/hashes.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/mpi/allreduce.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/utils/block_utils.h>

namespace hydra::indexing {

template <typename bit_t>
SpinhalfMPIIndexingSz<bit_t>::SpinhalfMPIIndexingSz(int n_sites, int n_up)
    : n_sites_(n_sites), n_up_(n_up), n_prefix_bits_(n_sites / 2),
      n_postfix_bits_(n_sites - n_prefix_bits_), empty_lintable_(), size_(0),
      size_transpose_(0) {
  using namespace combinatorics;

  utils::check_nup_spinhalf(n_sites, n_up, "SpinhalfMPIIndexingSz");

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

  // Determine the valid prefixes that belong to my process
  for (auto prefix : Subsets<bit_t>(n_prefix_bits_)) {
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
  size_transpose_ = 0;
  for (auto postfix : Subsets<bit_t>(n_postfix_bits_)) {
    int n_up_postfix = bitops::popcnt(postfix);
    int n_up_prefix = n_up_ - n_up_postfix;
    if ((n_up_prefix < 0) || (n_up_prefix > n_prefix_bits_))
      continue;

    // only keep "random" postfixes
    if (process(postfix) != mpi_rank_)
      continue;

    postfix_begin_[postfix] = size_transpose_;
    size_transpose_ += combinatorics::binomial(n_prefix_bits_, n_up_prefix);
    postfixes_.push_back(postfix);
  }

  size_max_ = std::max(size_, size_transpose_);

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

  // Compute total dimension and check whether correct
  dim_ = 0;
  mpi::Allreduce(&size_, &dim_, 1, MPI_SUM, MPI_COMM_WORLD);
  assert(dim_ == combinatorics::binomial(n_sites_, n_up_));
  mpi::Allreduce(&size_transpose_, &dim_, 1, MPI_SUM, MPI_COMM_WORLD);
  assert(dim_ == combinatorics::binomial(n_sites_, n_up_));
}

template <typename bit_t>
int SpinhalfMPIIndexingSz<bit_t>::process(bit_t prepostfix) const {
  return combinatorics::hash_fnv1(prepostfix) % mpi_size_;
}

template <typename bit_t>
gsl::span<bit_t const> SpinhalfMPIIndexingSz<bit_t>::prefixes() const {
  return {prefixes_.data(), prefixes_.size()};
}

template <typename bit_t>
gsl::span<bit_t const>
SpinhalfMPIIndexingSz<bit_t>::postfixes(bit_t prefix) const {
  int n_up_prefix = bitops::popcnt(prefix);
  int n_up_postfix = n_up_ - n_up_prefix;
  // assert(n_up_postfix >= 0);
  if ((n_up_postfix >= 0) && (n_up_postfix <= n_up_)) {
    return postfix_states_[n_up_postfix];
  } else {
    return gsl::span<bit_t const>();
  }
}

template <typename bit_t>
LinTable<bit_t> const &
SpinhalfMPIIndexingSz<bit_t>::postfix_indexing(bit_t prefix) const {
  int n_up_prefix = bitops::popcnt(prefix);
  int n_up_postfix = n_up_ - n_up_prefix;
  // lila::Log("nup: {} nuppre: {} nuppost: {}", n_up_, n_up_prefix,
  // n_up_postfix); assert(n_up_postfix >= 0);
  if ((n_up_postfix >= 0) && (n_up_postfix <= n_up_)) {
    return postfix_lintables_[n_up_postfix];
  } else {
    return empty_lintable_;
  }
}

template <typename bit_t>
gsl::span<bit_t const> SpinhalfMPIIndexingSz<bit_t>::postfixes() const {
  return {postfixes_.data(), postfixes_.size()};
}

template <typename bit_t>
gsl::span<bit_t const>
SpinhalfMPIIndexingSz<bit_t>::prefixes(bit_t postfix) const {
  int n_up_postfix = bitops::popcnt(postfix);
  int n_up_prefix = n_up_ - n_up_postfix;
  // assert(n_up_prefix >= 0);
  // return prefix_states_[n_up_prefix];
  if ((n_up_prefix >= 0) && (n_up_prefix <= n_up_)) {
    return prefix_states_[n_up_prefix];
  } else {
    return gsl::span<bit_t const>();
  }
}

template <typename bit_t>
LinTable<bit_t> const &
SpinhalfMPIIndexingSz<bit_t>::prefix_indexing(bit_t postfix) const {
  int n_up_postfix = bitops::popcnt(postfix);
  int n_up_prefix = n_up_ - n_up_postfix;
  // assert(n_up_prefix >= 0);
  // return prefix_lintables_[n_up_prefix];

  if ((n_up_prefix >= 0) && (n_up_prefix <= n_up_)) {
    return prefix_lintables_[n_up_prefix];
  } else {
    return empty_lintable_;
  }
}
template class SpinhalfMPIIndexingSz<uint16_t>;
template class SpinhalfMPIIndexingSz<uint32_t>;
template class SpinhalfMPIIndexingSz<uint64_t>;

} // namespace hydra::indexing
