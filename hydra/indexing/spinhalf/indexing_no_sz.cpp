#include "indexing_no_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::indexing::spinhalf {

template <typename bit_t>
IndexingNoSz<bit_t>::IndexingNoSz(int n_sites)
    : n_sites_(n_sites), size_(pow(2, n_sites)), begin_(0), end_(size_) {
  assert(n_sites_ >= 0);
}

#ifdef HYDRA_ENABLE_OPENMP
template <class bit_t>
typename IndexingNoSz<bit_t>::iterator_t
IndexingNoSz<bit_t>::thread_begin() const {
  idx_t start = omp::get_omp_start(size_);
  return iterator_t(start);
}

template <class bit_t>
typename IndexingNoSz<bit_t>::iterator_t
IndexingNoSz<bit_t>::thread_end() const {
  idx_t end = omp::get_omp_end(size_);
  return iterator_t(end);
}
#endif

template class IndexingNoSz<uint16_t>;
template class IndexingNoSz<uint32_t>;
template class IndexingNoSz<uint64_t>;

} // namespace hydra::indexing::spinhalf
