#include "combinations_indexing.h"

#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing {

template <class bit_t>
CombinationsIndexing<bit_t>::CombinationsIndexing(int n, int k)
  : n_(n), k_(k), size_(combinatorics::binomial(n, k)), lin_table_(n, k) {
  if (k > n) {
    lila::Log.err("Error constructing CombinationsIndexing: k>n");
  } else if (k < 0) {
    lila::Log.err("Error constructing CombinationsIndexing: k<0");
  } else if (n < 0) {
    lila::Log.err("Error constructing CombinationsIndexing: n<0");
  }
}

template class CombinationsIndexing<uint16_t>;
template class CombinationsIndexing<uint32_t>;
template class CombinationsIndexing<uint64_t>;
} // namespace hydra::indexing
