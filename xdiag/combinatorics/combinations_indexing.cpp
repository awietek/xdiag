#include "combinations_indexing.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
CombinationsIndexing<bit_t>::CombinationsIndexing(int64_t n, int64_t k)
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)), lin_table_(n, k) {
  if (k > n) {
    Log.err("Error constructing CombinationsIndexing: k>n");
  } else if (k < 0) {
    Log.err("Error constructing CombinationsIndexing: k<0");
  } else if (n < 0) {
    Log.err("Error constructing CombinationsIndexing: n<0");
  }
}

template class CombinationsIndexing<uint16_t>;
template class CombinationsIndexing<uint32_t>;
template class CombinationsIndexing<uint64_t>;
} // namespace xdiag::combinatorics
