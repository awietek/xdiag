#include "indexing_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations_index.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::indexing::spinhalf {

template <typename bit_t>
IndexingSz<bit_t>::IndexingSz(int n_sites, int nup)
    : n_sites_(n_sites), n_up_(nup), lintable_(n_sites, nup),
      states_(combinatorics::binomial(n_sites, nup)), size_(states_.size()),
      begin_(n_sites, nup, 0), end_(n_sites, nup, size_) {
  utils::check_nup_spinhalf(n_sites, nup, "Spinhalf");

#ifdef _OPENMP
  for (auto [state, idx] :
       combinatorics::CombinationsIndexThread<bit_t>(n_sites, nup)) {
    states_[idx] = state;
  }
#else
  for (auto [state, idx] :
       combinatorics::CombinationsIndex<bit_t>(n_sites, nup)) {
    states_[idx] = state;
  }
#endif
}

template <typename bit_t>
bool IndexingSz<bit_t>::operator==(IndexingSz<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_);
}

template <typename bit_t>
bool IndexingSz<bit_t>::operator!=(IndexingSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class IndexingSz<uint16_t>;
template class IndexingSz<uint32_t>;
template class IndexingSz<uint64_t>;

} // namespace hydra::indexing::spinhalf
