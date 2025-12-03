// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "subsets_table.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <typename bit_t>
SubsetsTable<bit_t>::SubsetsTable(int64_t n) try
    : n_(n), size_((int64_t)1 << n) {
  if (n < 0) {
    XDIAG_THROW("Error constructing SubsetsTable: n<0");
  }
}
XDIAG_CATCH

template <typename bit_t> int64_t SubsetsTable<bit_t>::size() const {
  return size_;
}

template <typename bit_t> Subsets<bit_t> SubsetsTable<bit_t>::states() const {
  return Subsets<bit_t>(n_);
}

template class SubsetsTable<uint16_t>;
template class SubsetsTable<uint32_t>;
template class SubsetsTable<uint64_t>;

} // namespace xdiag::combinatorics
