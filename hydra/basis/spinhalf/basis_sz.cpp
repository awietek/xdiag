#include "basis_sz.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations_index.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::basis::spinhalf {

template <typename bit_t>
BasisSz<bit_t>::BasisSz(int64_t n_sites, int64_t nup)
    : n_sites_(n_sites), n_up_(nup), lintable_(n_sites, nup),
      states_(combinatorics::binomial(n_sites, nup)), size_(states_.size()),
      begin_(n_sites, nup, 0), end_(n_sites, nup, size_) {
  if ((nup < 0) || (nup > n_sites)) {
    throw(std::invalid_argument("Invalid value of nup"));
  } else if (n_sites < 0) {
    throw(std::invalid_argument("n_sites < 0"));
  }

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
combinatorics::CombinationsIterator<bit_t> BasisSz<bit_t>::begin() const {
  return begin_;
}
template <typename bit_t>
combinatorics::CombinationsIterator<bit_t> BasisSz<bit_t>::end() const {
  return end_;
}

template <typename bit_t> int64_t BasisSz<bit_t>::size() const { return size_; }
template <typename bit_t> int64_t BasisSz<bit_t>::n_sites() const {
  return n_sites_;
}
template <typename bit_t> int64_t BasisSz<bit_t>::n_up() const { return n_up_; }

template <typename bit_t>
bool BasisSz<bit_t>::operator==(BasisSz<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_up_ == rhs.n_up_);
}

template <typename bit_t>
bool BasisSz<bit_t>::operator!=(BasisSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BasisSz<uint16_t>;
template class BasisSz<uint32_t>;
template class BasisSz<uint64_t>;

} // namespace hydra::basis::spinhalf
