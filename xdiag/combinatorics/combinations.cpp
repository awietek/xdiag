#include "combinations.h"

#include <xdiag/combinatorics/bit_patterns.h>
#include <xdiag/utils/logger.h>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.h>
#endif

namespace xdiag::combinatorics {

template <class bit_t>
Combinations<bit_t>::Combinations(int64_t n, int64_t k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    Log.err("Error constructing Combinations: k>n");
  } else if (k < 0) {
    Log.err("Error constructing Combinations: k<0");
  } else if (n < 0) {
    Log.err("Error constructing Combinations: n<0");
  } else {
    begin_ = CombinationsIterator<bit_t>(n, k, (int64_t)0);
    end_ = CombinationsIterator<bit_t>(n, k, (int64_t)binomial(n, k));
  }
}

template class Combinations<uint16_t>;
template class Combinations<uint32_t>;
template class Combinations<uint64_t>;

template <class bit_t>
CombinationsIterator<bit_t>::CombinationsIterator(int64_t n, int64_t k,
                                                  int64_t idx)
    : current_(get_nth_pattern<bit_t>(idx, n, k)), idx_(idx) {}

template class CombinationsIterator<uint16_t>;
template class CombinationsIterator<uint32_t>;
template class CombinationsIterator<uint64_t>;

#ifdef _OPENMP
template <class bit_t>
CombinationsThread<bit_t>::CombinationsThread(int64_t n, int64_t k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    Log.err("Error constructing Combinations: k>n");
  } else if (k < 0) {
    Log.err("Error constructing Combinations: k<0");
  } else if (n < 0) {
    Log.err("Error constructing Combinations: n<0");
  } else {
    auto [begin_idx, end_idx] = omp::get_omp_combinations_start_end(n, k);
    begin_ = CombinationsIterator<bit_t>(n, k, begin_idx);
    end_ = CombinationsIterator<bit_t>(n, k, end_idx);
  }
}

template class CombinationsThread<uint16_t>;
template class CombinationsThread<uint32_t>;
template class CombinationsThread<uint64_t>;
#endif

} // namespace xdiag::combinatorics
