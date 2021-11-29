#include "combinations.h"

#include <omp.h>

#include <lila/utils/logger.h>

#include <hydra/combinatorics/bit_patterns.h>

namespace hydra::combinatorics {

template <class bit_t>
Combinations<bit_t>::Combinations(int n, int k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    lila::Log.err("Error constructing Combinations: k>n");
  } else if (k < 0) {
    lila::Log.err("Error constructing Combinations: k<0");
  } else if (n < 0) {
    lila::Log.err("Error constructing Combinations: n<0");
  } else {
    bit_t state = ((bit_t)1 << k) - 1;
    begin_ = CombinationsIterator<bit_t>(state, (idx_t)0);
    end_ = CombinationsIterator<bit_t>(state,
                                       (idx_t)combinatorics::binomial(n, k));
  }
}

template <class bit_t>
CombinationsIterator<bit_t>::CombinationsIterator(bit_t state, idx_t idx)
    : current_(state), idx_(idx) {}

template class Combinations<uint16_t>;
template class Combinations<uint32_t>;
template class Combinations<uint64_t>;

template class CombinationsIterator<uint16_t>;
template class CombinationsIterator<uint32_t>;
template class CombinationsIterator<uint64_t>;

} // namespace hydra::combinatorics
