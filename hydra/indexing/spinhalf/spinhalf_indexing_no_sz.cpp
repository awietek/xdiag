#include "spinhalf_indexing_no_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing {

template <typename bit_t>
SpinhalfIndexingNoSz<bit_t>::SpinhalfIndexingNoSz(int n_sites)
    : n_sites_(n_sites), size_(pow(2, n_sites)) {
  assert(n_sites_ >= 0);
}

template class SpinhalfIndexingNoSz<uint16_t>;
template class SpinhalfIndexingNoSz<uint32_t>;
template class SpinhalfIndexingNoSz<uint64_t>;

} // namespace hydra::indexing
