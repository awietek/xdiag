#include "combinations_index.h"

#include <xdiag/combinatorics/bit_patterns.h>
#include <xdiag/utils/logger.h>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.h>
#endif

namespace xdiag::combinatorics {

template <class bit_t>
CombinationsIndex<bit_t>::CombinationsIndex(int64_t n, int64_t k)
    : n_(n), k_(k), size_(binomial(n, k)) {
  if (k > n) {
    Log.err("Error constructing CombinationsIndex: k>n");
  } else if (k < 0) {
    Log.err("Error constructing CombinationsIndex: k<0");
  } else if (n < 0) {
    Log.err("Error constructing CombinationsIndex: n<0");
  } else {
    begin_ = CombinationsIndexIterator<bit_t>(n, k, 0);
    end_ = CombinationsIndexIterator<bit_t>(n, k, binomial(n, k));
  }
}

template <class bit_t>
CombinationsIndexIterator<bit_t>::CombinationsIndexIterator(int64_t n,
                                                            int64_t k,
                                                            int64_t idx)
    : current_(get_nth_pattern<bit_t>(idx, n, k)), idx_(idx) {}

template class CombinationsIndex<uint16_t>;
template class CombinationsIndex<uint32_t>;
template class CombinationsIndex<uint64_t>;

template class CombinationsIndexIterator<uint16_t>;
template class CombinationsIndexIterator<uint32_t>;
template class CombinationsIndexIterator<uint64_t>;

#ifdef _OPENMP
template <class bit_t>
CombinationsIndexThread<bit_t>::CombinationsIndexThread(int64_t n, int64_t k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    Log.err("Error constructing CombinationsIndexThread: k>n");
  } else if (k < 0) {
    Log.err("Error constructing CombinationsIndexThread: k<0");
  } else if (n < 0) {
    Log.err("Error constructing CombinationsIndexThread: n<0");
  } else {
    auto [begin_idx, end_idx] = omp::get_omp_combinations_start_end(n, k);
    begin_ = CombinationsIndexIterator<bit_t>(n, k, begin_idx);
    end_ = CombinationsIndexIterator<bit_t>(n, k, end_idx);
  }
}

template class CombinationsIndexThread<uint16_t>;
template class CombinationsIndexThread<uint32_t>;
template class CombinationsIndexThread<uint64_t>;
#endif

} // namespace xdiag::combinatorics
