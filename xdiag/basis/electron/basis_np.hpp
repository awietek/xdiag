// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::electron {

template <typename bit_tt> class BasisNpIterator;

using namespace combinatorics;

template <typename bit_tt> class BasisNp {
public:
  using bit_t = bit_tt;
  using iterator_t = BasisNpIterator<bit_t>;

  BasisNp() = default;
  BasisNp(int nsites, int nup, int ndn);

  int64_t size() const;
  int64_t dim() const;
  iterator_t begin() const;
  iterator_t end() const;

  inline int64_t index(bit_t ups, bit_t dns) const {
    return index_ups(ups) * size_dns_ + index_dns(dns);
  }

  int nsites() const;
  int nup() const;
  int ndn() const;
  static constexpr bool np_conserved() { return true; }

private:
  int nsites_;
  int nup_;
  int ndn_;

  int64_t size_ups_;
  int64_t size_dns_;
  int64_t size_;

  combinatorics::LinTable<bit_t> lintable_ups_;
  combinatorics::LinTable<bit_t> lintable_dns_;

public:
  Combinations<bit_t> states_ups() const;
  Combinations<bit_t> states_dns() const;
  CombinationsIndex<bit_t> states_indices_ups() const;
  CombinationsIndex<bit_t> states_indices_dns() const;

#ifdef _OPENMP
  CombinationsThread<bit_t> states_ups_thread() const;
  CombinationsThread<bit_t> states_dns_thread() const;
  CombinationsIndexThread<bit_t> states_indices_ups_thread() const;
  CombinationsIndexThread<bit_t> states_indices_dns_thread() const;
#endif

  int64_t size_ups() const;
  int64_t size_dns() const;
  inline int64_t index_ups(bit_t ups) const { return lintable_ups_.index(ups); }
  inline int64_t index_dns(bit_t dns) const { return lintable_dns_.index(dns); }
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
  bit_t begindns_;
  bit_t end_dns_;
  bit_t ups_;
  bit_t dns_;
};
} // namespace xdiag::basis::electron
