#include "combinations_index.h"

#include <omp.h>

#include <lila/utils/logger.h>

#include <hydra/combinatorics/bit_patterns.h>

namespace hydra::combinatorics {

template <class bit_t>
CombinationsIndex<bit_t>::CombinationsIndex(int n, int k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    lila::Log.err("Error constructing CombinationsIndex: k>n");
  } else if (k < 0) {
    lila::Log.err("Error constructing CombinationsIndex: k<0");
  } else if (n < 0) {
    lila::Log.err("Error constructing CombinationsIndex: n<0");
  } else {
    bit_t state = ((bit_t)1 << k) - 1;
    begin_ = CombinationsIndexIterator<bit_t>(state, (idx_t)0);
    end_ = CombinationsIndexIterator<bit_t>(
        state, (idx_t)combinatorics::binomial(n, k));
  }
}

template <class bit_t>
CombinationsIndexIterator<bit_t>::CombinationsIndexIterator(bit_t state,
                                                            idx_t idx)
    : current_(state), idx_(idx) {}

template class CombinationsIndex<uint16_t>;
template class CombinationsIndex<uint32_t>;
template class CombinationsIndex<uint64_t>;

template class CombinationsIndexIterator<uint16_t>;
template class CombinationsIndexIterator<uint32_t>;
template class CombinationsIndexIterator<uint64_t>;

template <typename bit_t>
std::pair<CombinationsIndexIterator<bit_t>, CombinationsIndexIterator<bit_t>>
get_omp_combinations_start_end(int n, int k) {

#ifndef _OPENMP
  bit_t state = ((bit_t)1 << k) - 1;
  auto start = CombinationsIndexIterator<bit_t>(state, (idx_t)0);
  auto end = CombinationsIndexIterator<bit_t>(
      state, (idx_t)combinatorics::binomial(n, k));
#else
  idx_t myid = omp_get_thread_num();
  idx_t rank = omp_get_num_threads();
  idx_t size = binomial(n, k);
  idx_t chunksize = size / rank;

  idx_t start_idx = myid * chunksize;
  idx_t end_idx = (myid == rank - 1) ? size : (myid + 1) * chunksize;

  bit_t start_state = get_nth_pattern<bit_t>(start_idx, n, k);
  bit_t end_state = get_nth_pattern<bit_t>(end_idx, n, k);

  auto start = CombinationsIndexIterator<bit_t>(start_state, start_idx);
  auto end = CombinationsIndexIterator<bit_t>(end_state, end_idx);
#endif
  return {start, end};
}

template std::pair<CombinationsIndexIterator<uint16_t>,
                   CombinationsIndexIterator<uint16_t>>
get_omp_combinations_start_end<uint16_t>(int n, int k);
template std::pair<CombinationsIndexIterator<uint32_t>,
                   CombinationsIndexIterator<uint32_t>>
get_omp_combinations_start_end<uint32_t>(int n, int k);
template std::pair<CombinationsIndexIterator<uint64_t>,
                   CombinationsIndexIterator<uint64_t>>
get_omp_combinations_start_end<uint64_t>(int n, int k);

template <class bit_t>
CombinationsIndexThread<bit_t>::CombinationsIndexThread(int n, int k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    lila::Log.err("Error constructing CombinationsIndexThread: k>n");
  } else if (k < 0) {
    lila::Log.err("Error constructing CombinationsIndexThread: k<0");
  } else if (n < 0) {
    lila::Log.err("Error constructing CombinationsIndexThread: n<0");
  } else {
    std::tie(begin_, end_) = get_omp_combinations_start_end<bit_t>(n, k);
  }
}

template class CombinationsIndexThread<uint16_t>;
template class CombinationsIndexThread<uint32_t>;
template class CombinationsIndexThread<uint64_t>;

} // namespace hydra::combinatorics
