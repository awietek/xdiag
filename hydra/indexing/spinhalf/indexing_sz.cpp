#include "indexing_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::indexing::spinhalf {

template <typename bit_t>
IndexingSz<bit_t>::IndexingSz(int n_sites, int nup)
    : n_sites_(n_sites), n_up_(nup), lintable_(n_sites, nup),
      size_(combinatorics::binomial(n_sites, nup)), begin_(n_sites, nup, 0),
      end_(n_sites, nup, size_) {
  utils::check_nup_spinhalf(n_sites, nup, "Spinhalf");
}

#ifdef HYDRA_ENABLE_OPENMP
template <class bit_t>
typename IndexingSz<bit_t>::iterator_t IndexingSz<bit_t>::thread_begin() const {
  idx_t start = omp::get_omp_start(size_);
  return iterator_t(n_sites_, n_up_, start);
}

template <class bit_t>
typename IndexingSz<bit_t>::iterator_t IndexingSz<bit_t>::thread_end() const {
  idx_t end = omp::get_omp_end(size_);
  return iterator_t(n_sites_, n_up_, end);
}
#endif

template class IndexingSz<uint16_t>;
template class IndexingSz<uint32_t>;
template class IndexingSz<uint64_t>;

} // namespace hydra::indexing::spinhalf
