// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "combinations_table.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
CombinationsTable<bit_t>::CombinationsTable(int64_t n, int64_t k) try
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)), lin_table_(n, k) {
  if (k > n) {
    XDIAG_THROW("Error constructing CombinationsTable: k>n");
  } else if (k < 0) {
    XDIAG_THROW("Error constructing CombinationsTable: k<0");
  } else if (n < 0) {
    XDIAG_THROW("Error constructing CombinationsTable: n<0");
  }
}
XDIAG_CATCH

template <class bit_t>
bit_t CombinationsTable<bit_t>::operator[](int64_t idx) const {
  return combinatorics::get_nth_pattern<bit_t>(idx, n_, k_);
}

template <class bit_t> int64_t CombinationsTable<bit_t>::size() const {
  return size_;
}

template <class bit_t>
Combinations<bit_t> CombinationsTable<bit_t>::states() const {
  return Combinations<bit_t>(n_, k_);
}

template <class bit_t>
bool CombinationsTable<bit_t>::operator==(
    CombinationsTable<bit_t> const &rhs) const {
  return (n_ == rhs.n_) && (k_ == rhs.k_);
}

template <class bit_t>
bool CombinationsTable<bit_t>::operator!=(
    CombinationsTable<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class CombinationsTable<uint16_t>;
template class CombinationsTable<uint32_t>;
template class CombinationsTable<uint64_t>;

template <class bit_t>
CombinationsTable<bit_t> table(Combinations<bit_t> const &states) {
  return CombinationsTable<bit_t>(states.n(), states.k());
}

template CombinationsTable<uint16_t> table(Combinations<uint16_t> const &);
template CombinationsTable<uint32_t> table(Combinations<uint32_t> const &);
template CombinationsTable<uint64_t> table(Combinations<uint64_t> const &);
} // namespace xdiag::combinatorics
