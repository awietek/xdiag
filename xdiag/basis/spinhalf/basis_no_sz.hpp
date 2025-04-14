// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/subsets_index.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_tt> class BasisNoSz {
public:
  using bit_t = bit_tt;
  using iterator_t = combinatorics::SubsetsIterator<bit_t>;

  BasisNoSz() = default;
  BasisNoSz(int64_t nsites);

  iterator_t begin() const;
  iterator_t end() const;
  int64_t dim() const;
  int64_t size() const;

  inline int64_t index(bit_t spins) const { return (int64_t)spins; }
  inline bit_t state(int64_t index) const { return (bit_t)index; }

  int64_t nsites() const;
  static constexpr bool sz_conserved() { return false; }

  bool operator==(BasisNoSz const &rhs) const;
  bool operator!=(BasisNoSz const &rhs) const;

private:
  int64_t nsites_;
  int64_t size_;
  iterator_t begin_, end_;
};

} // namespace xdiag::basis::spinhalf
