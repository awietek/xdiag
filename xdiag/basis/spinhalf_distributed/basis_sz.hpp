// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/extern/gsl/span>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>
#include <xdiag/parallel/mpi/comm_pattern.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>
#include <xdiag/random/hash_functions.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename bit_tt> class BasisSzIterator;

template <typename bit_tt> class BasisSz {
public:
  using bit_t = bit_tt;
  using iterator_t = BasisSzIterator<bit_t>;

  BasisSz() = default;
  BasisSz(int64_t nsites, int64_t nup);

  int64_t nsites() const;
  int64_t nup() const;

  int64_t n_prefix_bits() const;
  int64_t n_postfix_bits() const;

  int64_t dim() const;
  int64_t size() const;
  int64_t size_transpose() const;
  int64_t size_max() const;
  int64_t size_min() const;
  iterator_t begin() const;
  iterator_t end() const;
  int64_t index(bit_t spins) const;

  std::vector<bit_t> const &prefixes() const;
  int64_t prefix_begin(bit_t prefix) const;
  combinatorics::LinTable<bit_t> const &postfix_lintable(bit_t prefix) const;
  std::vector<bit_t> const &postfix_states(bit_t prefix) const;

  std::vector<bit_t> const &postfixes() const;
  int64_t postfix_begin(bit_t postfix) const;
  combinatorics::LinTable<bit_t> const &prefix_lintable(bit_t postfix) const;
  std::vector<bit_t> const &prefix_states(bit_t postfix) const;

  inline int rank(bit_t spins) const {
    return (int)(random::hash_div3(spins) % mpi_size_);
  };

  mpi::CommPattern &comm_pattern() const;
  mpi::Communicator transpose_communicator(bool reverse) const;

  bool operator==(BasisSz const &rhs) const;
  bool operator!=(BasisSz const &rhs) const;

private:
  int64_t nsites_;
  int64_t nup_;

  int n_prefix_bits_;
  int n_postfix_bits_;

  int64_t dim_;
  int64_t size_;
  int64_t size_transpose_;
  int64_t size_max_;
  int64_t size_min_;

  int mpi_rank_;
  int mpi_size_;

  std::vector<bit_t> prefixes_;
  std::unordered_map<bit_t, int64_t> prefix_begin_;
  std::vector<combinatorics::LinTable<bit_t>> postfix_lintables_;
  std::vector<std::vector<bit_t>> postfix_states_;

  std::vector<bit_t> postfixes_;
  std::unordered_map<bit_t, int64_t> postfix_begin_;
  std::vector<combinatorics::LinTable<bit_t>> prefix_lintables_;
  std::vector<std::vector<bit_t>> prefix_states_;

  mutable mpi::CommPattern comm_pattern_;
  mpi::Communicator transpose_communicator_;
  mpi::Communicator transpose_communicator_reverse_;
};

template <typename bit_tt> class BasisSzIterator {
public:
  using bit_t = bit_tt;
  BasisSzIterator() = default;
  BasisSzIterator(BasisSz<bit_t> const& basis, bool begin);
  BasisSzIterator<bit_t> &operator++();
  bit_t operator*() const;
  bool operator!=(BasisSzIterator<bit_t> const &rhs) const;

private:
  BasisSz<bit_t> const &basis_;
  bit_t prefix_;
  int64_t prefix_idx_;
  int64_t postfix_idx_;
};

} // namespace xdiag::basis::spinhalf_distributed
#endif
