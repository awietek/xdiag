#include "basis_sz.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::basis::spinhalf {

template <typename bit_t>
BasisSz<bit_t>::BasisSz(int64_t n_sites, int64_t nup) try
    : n_sites_(n_sites), n_up_(nup), lintable_(n_sites, nup),
      states_(combinatorics::binomial(n_sites, nup)), size_(states_.size()),
      begin_(n_sites, nup, 0), end_(n_sites, nup, size_) {
  check_n_sites_work_with_bits<bit_t>(n_sites);

  if ((nup < 0) || (nup > n_sites)) {
    XDIAG_THROW("Invalid value of nup");
  } else if (n_sites < 0) {
    XDIAG_THROW("n_sites < 0");
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
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t>
combinatorics::CombinationsIterator<bit_t> BasisSz<bit_t>::begin() const {
  return begin_;
}
template <typename bit_t>
combinatorics::CombinationsIterator<bit_t> BasisSz<bit_t>::end() const {
  return end_;
}

template <typename bit_t> int64_t BasisSz<bit_t>::dim() const { return size_; }
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

template class BasisSz<uint32_t>;
template class BasisSz<uint64_t>;

} // namespace xdiag::basis::spinhalf
