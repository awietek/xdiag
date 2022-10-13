#include "indexing_no_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra::indexing::spinhalf {

template <typename bit_t>
IndexingNoSz<bit_t>::IndexingNoSz(int n_sites)
    : n_sites_(n_sites), size_(pow(2, n_sites)), begin_(0), end_(size_) {
  assert(n_sites_ >= 0);
}

template <typename bit_t>
bool IndexingNoSz<bit_t>::operator==(IndexingNoSz<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_);
}

template <typename bit_t>
bool IndexingNoSz<bit_t>::operator!=(IndexingNoSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class IndexingNoSz<uint16_t>;
template class IndexingNoSz<uint32_t>;
template class IndexingNoSz<uint64_t>;

} // namespace hydra::indexing::spinhalf
