#include "basis_sz.hpp"

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::basis::spinhalf {

template <typename bit_t>
BasisSz<bit_t>::BasisSz(int64_t nsites, int64_t nup) try
    : nsites_(nsites), nup_(nup), lintable_(nsites, nup),
      states_(combinatorics::binomial(nsites, nup)), size_(states_.size()),
      begin_(nsites, nup, 0), end_(nsites, nup, size_) {
  check_nsites_work_with_bits<bit_t>(nsites);

  if ((nup < 0) || (nup > nsites)) {
    XDIAG_THROW("Invalid value of nup");
  } else if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  }

#ifdef _OPENMP
  for (auto [state, idx] :
       combinatorics::CombinationsIndexThread<bit_t>(nsites, nup)) {
    states_[idx] = state;
  }
#else
  for (auto [state, idx] :
       combinatorics::CombinationsIndex<bit_t>(nsites, nup)) {
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
template <typename bit_t> int64_t BasisSz<bit_t>::nsites() const {
  return nsites_;
}
template <typename bit_t> int64_t BasisSz<bit_t>::nup() const { return nup_; }

template <typename bit_t>
bool BasisSz<bit_t>::operator==(BasisSz<bit_t> const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_);
}

template <typename bit_t>
bool BasisSz<bit_t>::operator!=(BasisSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BasisSz<uint32_t>;
template class BasisSz<uint64_t>;

} // namespace xdiag::basis::spinhalf
