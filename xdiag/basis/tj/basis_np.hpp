// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::tj {

template <typename bit_tt> class BasisNpIterator;

template <typename bit_tt> class BasisNp {
public:
  using bit_t = bit_tt;
  using iterator_t = BasisNpIterator<bit_t>;

  BasisNp() = default;
  BasisNp(int64_t nsites, int64_t nup, int64_t ndn);

  int64_t dim() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(BasisNp const &rhs) const;
  bool operator!=(BasisNp const &rhs) const;
  
  inline int64_t index(bit_t ups, bit_t dns) const {
    bit_t dncs = bits::extract<bit_t>(dns, (~ups) & sitesmask_);
    return index_ups(ups) * size_dncs_ + index_dncs(dncs);
  }

  int64_t nsites() const;
  int64_t nup() const;
  int64_t ndn() const;
  static constexpr bool np_conserved() { return true; }

private:
  int64_t nsites_;
  int64_t nup_;
  int64_t ndn_;

  int64_t size_ups_;
  int64_t size_dncs_;
  int64_t size_;
  bit_t sitesmask_;

  combinatorics::LinTable<bit_t> lintable_ups_;
  combinatorics::LinTable<bit_t> lintable_dncs_;

public:
  int64_t size_ups() const;
  int64_t size_dncs() const;
  int64_t ups_offset(int64_t idx_ups) const;

  combinatorics::Combinations<bit_t> states_ups() const;
  combinatorics::Combinations<bit_t> states_dncs(bit_t ups) const;
  combinatorics::CombinationsIndex<bit_t> states_indices_ups() const;
  combinatorics::CombinationsIndex<bit_t> states_indices_dncs(bit_t ups) const;

#ifdef _OPENMP
  combinatorics::CombinationsThread<bit_t> states_ups_thread() const;
  combinatorics::CombinationsThread<bit_t> states_dncs_thread(bit_t ups) const;
  combinatorics::CombinationsIndexThread<bit_t>
  states_indices_ups_thread() const;
  combinatorics::CombinationsIndexThread<bit_t>
  states_indices_dncs_thread(bit_t ups) const;
#endif

  inline int64_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }
  inline int64_t index_dncs(bit_t dncs) const {
    return lintable_dncs_.index(dncs);
  }
};

template <typename bit_tt> class BasisNpIterator {
public:
  using bit_t = bit_tt;
  BasisNpIterator() = default;
  BasisNpIterator(int64_t nsites, int64_t nup, int64_t ndn, bool begin);
  BasisNpIterator &operator++();
  std::pair<bit_t, bit_t> operator*() const;
  bool operator!=(BasisNpIterator<bit_t> const &rhs) const;

private:
  bit_t sitesmask_;
  bit_t begindns_;
  bit_t end_dns_;
  bit_t ups_;
  bit_t dns_;
};

} // namespace xdiag::basis::tj
