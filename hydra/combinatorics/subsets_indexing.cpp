#include "subsets_indexing.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/utils/logger.h>

namespace hydra::combinatorics {

template <class bit_t>
SubsetsIndexing<bit_t>::SubsetsIndexing(int n) : n_(n), size_((idx_t)1 << n) {
  if (n < 0) {
    Log.err("Error constructing SubsetsIndexing: n<0");
  }
}

template class SubsetsIndexing<uint16_t>;
template class SubsetsIndexing<uint32_t>;
template class SubsetsIndexing<uint64_t>;

} // namespace hydra::combinatorics
