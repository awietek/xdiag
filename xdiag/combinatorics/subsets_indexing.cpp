// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "subsets_indexing.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
SubsetsIndexing<bit_t>::SubsetsIndexing(int64_t n)
    : n_(n), size_((int64_t)1 << n) {
  if (n < 0) {
    Log.err("Error constructing SubsetsIndexing: n<0");
  }
}

template class SubsetsIndexing<uint16_t>;
template class SubsetsIndexing<uint32_t>;
template class SubsetsIndexing<uint64_t>;

} // namespace xdiag::combinatorics
