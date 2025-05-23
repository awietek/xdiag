// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "subsets_index.hpp"

#include <xdiag/utils/logger.hpp>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::combinatorics {

template <class bit_t>
SubsetsIndex<bit_t>::SubsetsIndex(int64_t n) : n_(n), size_((int64_t)1 << n) {
  if (n < 0) {
    Log.err("Error constructing SubsetsIndex: n<0");
  } else {
    begin_ = SubsetsIndexIterator<bit_t>(0);
    end_ = SubsetsIndexIterator<bit_t>((int64_t)1 << n);
  }
}

template class SubsetsIndex<uint16_t>;
template class SubsetsIndex<uint32_t>;
template class SubsetsIndex<uint64_t>;

template <class bit_t>
SubsetsIndexIterator<bit_t>::SubsetsIndexIterator(int64_t idx)
    : current_((bit_t)idx) {}

template class SubsetsIndexIterator<uint16_t>;
template class SubsetsIndexIterator<uint32_t>;
template class SubsetsIndexIterator<uint64_t>;

#ifdef _OPENMP
template <class bit_t>
SubsetsIndexThread<bit_t>::SubsetsIndexThread(int64_t n)
    : n_(n), size_((int64_t)1 << n) {
  if (n < 0) {
    Log.err("Error constructing SubsetsIndex: n<0");
  } else {
    auto [begin_idx, end_idx] = omp::get_omp_subsets_start_end(n);
    begin_ = SubsetsIndexIterator<bit_t>(begin_idx);
    end_ = SubsetsIndexIterator<bit_t>(end_idx);
  }
}

template class SubsetsIndexThread<uint16_t>;
template class SubsetsIndexThread<uint32_t>;
template class SubsetsIndexThread<uint64_t>;
#endif

} // namespace xdiag::combinatorics
