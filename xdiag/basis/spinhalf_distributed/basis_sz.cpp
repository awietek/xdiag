// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_sz.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/parallel/mpi/allreduce.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename bit_t, class process_f>
static int64_t
fill_tables(int64_t nsites, int64_t nup, int64_t n_prefix_bits, process_f rank,
            std::vector<bit_t> &prefixes,
            std::unordered_map<bit_t, int64_t> &prefix_begin,
            std::vector<combinatorics::LinTable<bit_t>> &postfix_lintables,
            std::vector<std::vector<bit_t>> &postfix_states) {

  using combinatorics::binomial;
  using combinatorics::LinTable;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int64_t n_postfix_bits = nsites - n_prefix_bits;

  // Determine the valid prefixes that belong to my process
  int64_t size = 0;
  for (bit_t prefix : combinatorics::Subsets<bit_t>(n_prefix_bits)) {
    int nup_prefix = bits::popcnt(prefix);
    int nup_postfix = nup - nup_prefix;

    // Ignore impossible prefix configurations
    if ((nup_postfix < 0) || (nup_postfix > n_postfix_bits)) {
      continue;
    }

    // only keep "random" prefixes
    if (rank(prefix) != mpi_rank) {
      continue;
    }

    // Register a prefix
    prefix_begin[prefix] = size;
    size += binomial(n_postfix_bits, nup_postfix);
    prefixes.push_back(prefix);
  }

  // Create the lintables for postfix lookup
  postfix_states.resize(n_postfix_bits + 1);
  for (int nup_postfix = 0; nup_postfix <= n_postfix_bits; ++nup_postfix) {
    postfix_lintables.push_back(LinTable<bit_t>(n_postfix_bits, nup_postfix));

    // Precompute posfix configurations for given nup_postfix
    std::vector<bit_t> postfixes_nup(binomial(n_postfix_bits, nup_postfix), 0);
    int64_t pf_idx = 0;
    for (bit_t pf :
         combinatorics::Combinations<bit_t>(n_postfix_bits, nup_postfix)) {
      postfixes_nup[pf_idx++] = pf;
    }
    postfix_states[nup_postfix] = postfixes_nup;
  }
  return size;
}

template <typename bit_t>
BasisSz<bit_t>::BasisSz(int64_t nsites, int64_t nup)
    : nsites_(nsites), nup_(nup), n_prefix_bits_(nsites / 2),
      n_postfix_bits_(nsites - n_prefix_bits_) {
  using namespace combinatorics;
  check_nsites_work_with_bits<bit_t>(nsites_);

  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if (nup < 0) {
    XDIAG_THROW("nup < 0");
  } else if (nup > nsites) {
    XDIAG_THROW("nup > nsites");
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);

  dim_ = binomial(nsites, nup);
  size_ = fill_tables(
      nsites, nup, n_prefix_bits_, [this](bit_t spins) { return rank(spins); },
      prefixes_, prefix_begin_, postfix_lintables_, postfix_states_);

  size_transpose_ = fill_tables(
      nsites, nup, n_postfix_bits_, [this](bit_t spins) { return rank(spins); },
      postfixes_, postfix_begin_, prefix_lintables_, prefix_states_);

  // Compute max/min number of states stored locally
  int64_t size_max;
  int64_t size_max_transpose;
  mpi::Allreduce(&size_, &size_max, 1, MPI_MAX, MPI_COMM_WORLD);
  mpi::Allreduce(&size_transpose_, &size_max_transpose, 1, MPI_MAX,
                 MPI_COMM_WORLD);
  size_max_ = std::max(size_max, size_max_transpose);

  int64_t size_min;
  int64_t size_min_transpose;
  mpi::Allreduce(&size_, &size_min, 1, MPI_MIN, MPI_COMM_WORLD);
  mpi::Allreduce(&size_transpose_, &size_min_transpose, 1, MPI_MIN,
                 MPI_COMM_WORLD);
  size_min_ = std::min(size_min, size_min_transpose);

  // Check local sizes sum up to the actual dimension
  int64_t dim;
  mpi::Allreduce(&size_, &dim, 1, MPI_SUM, MPI_COMM_WORLD);
  assert(dim == dim_);

  int64_t dim_transpose;
  mpi::Allreduce(&size_transpose_, &dim_transpose, 1, MPI_SUM, MPI_COMM_WORLD);
  assert(dim_transpose == dim_);

  // Create the transpose communicator
  std::vector<int64_t> n_states_i_send(mpi_size_, 0);
  for (bit_t prefix : prefixes()) {
    for (bit_t postfix : postfix_states(prefix)) {
      int target_rank = rank(postfix);
      ++n_states_i_send[target_rank];
    }
  }
  transpose_communicator_ = mpi::Communicator(n_states_i_send);

  // Create the transpose communicator (reverse)
  std::vector<int64_t> n_states_i_send_reverse(mpi_size_, 0);
  for (bit_t postfix : postfixes()) {
    for (bit_t prefix : prefix_states(postfix)) {
      int target_rank = rank(prefix);
      ++n_states_i_send_reverse[target_rank];
    }
  }
  transpose_communicator_reverse_ = mpi::Communicator(n_states_i_send_reverse);
}

template <typename bit_t> int64_t BasisSz<bit_t>::nsites() const {
  return nsites_;
}
template <typename bit_t> int64_t BasisSz<bit_t>::nup() const { return nup_; }
template <typename bit_t> int64_t BasisSz<bit_t>::n_prefix_bits() const {
  return n_prefix_bits_;
};
template <typename bit_t> int64_t BasisSz<bit_t>::n_postfix_bits() const {
  return n_postfix_bits_;
}

template <typename bit_t> int64_t BasisSz<bit_t>::dim() const { return dim_; }
template <typename bit_t> int64_t BasisSz<bit_t>::size() const { return size_; }
template <typename bit_t> int64_t BasisSz<bit_t>::size_transpose() const {
  return size_transpose_;
}
template <typename bit_t> int64_t BasisSz<bit_t>::size_max() const {
  return size_max_;
}
template <typename bit_t> int64_t BasisSz<bit_t>::size_min() const {
  return size_min_;
}
template <typename bit_t>
typename BasisSz<bit_t>::iterator_t BasisSz<bit_t>::begin() const {
  return iterator_t(*this, true);
}
template <typename bit_t>
typename BasisSz<bit_t>::iterator_t BasisSz<bit_t>::end() const {
  return iterator_t(*this, false);
}
template <typename bit_t> int64_t BasisSz<bit_t>::index(bit_t spins) const {
  bit_t prefix = spins >> n_postfix_bits_;
  if (rank(prefix) != mpi_rank_) {
    return invalid_index;
  }
  int64_t offset = prefix_begin(prefix);
  auto const &lintable = postfix_lintable(prefix);
  bit_t postfix = spins & (((bit_t)1 << n_postfix_bits_) - 1);
  return offset + lintable.index(postfix);
}

template <typename bit_t>
std::vector<bit_t> const &BasisSz<bit_t>::prefixes() const {
  return prefixes_;
}
template <typename bit_t>
int64_t BasisSz<bit_t>::prefix_begin(bit_t prefix) const {
  auto it = prefix_begin_.find(prefix);
  if (it != prefix_begin_.end()) {
    return it->second;
  } else {
    return invalid_index;
  }
}
template <typename bit_t>
combinatorics::LinTable<bit_t> const &
BasisSz<bit_t>::postfix_lintable(bit_t prefix) const {
  int nup_prefix = bits::popcnt(prefix);
  int nup_postfix = nup_ - nup_prefix;
  return postfix_lintables_[nup_postfix];
}
template <typename bit_t>
std::vector<bit_t> const &BasisSz<bit_t>::postfix_states(bit_t prefix) const {
  int nup_prefix = bits::popcnt(prefix);
  int nup_postfix = nup_ - nup_prefix;
  return postfix_states_[nup_postfix];
}

template <typename bit_t>
std::vector<bit_t> const &BasisSz<bit_t>::postfixes() const {
  return postfixes_;
}
template <typename bit_t>
int64_t BasisSz<bit_t>::postfix_begin(bit_t postfix) const {
  auto it = postfix_begin_.find(postfix);
  if (it != postfix_begin_.end()) {
    return it->second;
  } else {
    return invalid_index;
  }
}
template <typename bit_t>
combinatorics::LinTable<bit_t> const &
BasisSz<bit_t>::prefix_lintable(bit_t postfix) const {
  int nup_postfix = bits::popcnt(postfix);
  int nup_prefix = nup_ - nup_postfix;
  return prefix_lintables_[nup_prefix];
}
template <typename bit_t>
std::vector<bit_t> const &BasisSz<bit_t>::prefix_states(bit_t postfix) const {
  int nup_postfix = bits::popcnt(postfix);
  int nup_prefix = nup_ - nup_postfix;
  return prefix_states_[nup_prefix];
}

template <typename bit_t>
mpi::CommPattern &BasisSz<bit_t>::comm_pattern() const {
  return comm_pattern_;
}

template <typename bit_t>
mpi::Communicator BasisSz<bit_t>::transpose_communicator(bool reverse) const {
  return reverse ? transpose_communicator_reverse_ : transpose_communicator_;
}

template class BasisSz<uint32_t>;
template class BasisSz<uint64_t>;

template <typename bit_t>
BasisSzIterator<bit_t>::BasisSzIterator(BasisSz<bit_t> const &basis, bool begin)
    : basis_(basis), prefix_idx_(begin ? 0 : basis.prefixes().size()),
      postfix_idx_(0) {
  if ((basis.prefixes().size() > 0) && begin) {
    prefix_ = basis.prefixes()[0];
  }
}

template <typename bit_t>
BasisSzIterator<bit_t> &BasisSzIterator<bit_t>::operator++() {
  ++postfix_idx_;
  if (postfix_idx_ == basis_.postfix_states(prefix_).size()) {
    postfix_idx_ = 0;
    ++prefix_idx_;
    if (prefix_idx_ != basis_.prefixes().size()) {
      prefix_ = basis_.prefixes()[prefix_idx_];
    }
  }
  return *this;
}

template <typename bit_t> bit_t BasisSzIterator<bit_t>::operator*() const {
  return (prefix_ << basis_.n_postfix_bits()) |
         basis_.postfix_states(prefix_)[postfix_idx_];
}

template <typename bit_t>
bool BasisSzIterator<bit_t>::operator!=(
    BasisSzIterator<bit_t> const &rhs) const {
  return (prefix_idx_ != rhs.prefix_idx_) || (postfix_idx_ != rhs.postfix_idx_);
}

template class BasisSzIterator<uint32_t>;
template class BasisSzIterator<uint64_t>;

} // namespace xdiag::basis::spinhalf_distributed
