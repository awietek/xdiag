// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <extern/gsl/span>

#include <xdiag/basis/basis.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/extract_deposit.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/mpi/comm_pattern.hpp>
#include <xdiag/mpi/communicator.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/utils/type_name.hpp>

namespace xdiag::basis {

template <typename bit_tt> class BasisSpinhalfDistributedIterator;

template <typename bit_tt>
class BasisSpinhalfDistributed
    : public BasisType<BasisSpinhalfDistributed<bit_tt>> {
public:
  using bit_t = bit_tt;
  using iterator_t = BasisSpinhalfDistributedIterator<bit_t>;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisSpinhalfDistributed<bit_t>>();

  BasisSpinhalfDistributed() = default;
  BasisSpinhalfDistributed(int64_t nsites, int64_t nup);

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

  bool operator==(BasisSpinhalfDistributed const &rhs) const;
  bool operator!=(BasisSpinhalfDistributed const &rhs) const;

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

template <typename bit_tt> class BasisSpinhalfDistributedIterator {
public:
  using bit_t = bit_tt;
  BasisSpinhalfDistributedIterator() = default;
  BasisSpinhalfDistributedIterator(BasisSpinhalfDistributed<bit_t> const &basis,
                                   bool begin);
  BasisSpinhalfDistributedIterator<bit_t> &operator++();
  bit_t operator*() const;
  bool operator!=(BasisSpinhalfDistributedIterator<bit_t> const &rhs) const;

private:
  BasisSpinhalfDistributed<bit_t> const &basis_;
  bit_t prefix_;
  int64_t prefix_idx_;
  int64_t postfix_idx_;
};

} // namespace xdiag::basis
