#include "basis_no_sz.hpp"

#include <xdiag/combinatorics/binomial.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_t>
BasisNoSz<bit_t>::BasisNoSz(int64_t n_sites)
    : n_sites_(n_sites), size_(pow(2, n_sites)), begin_(0), end_(size_) {
  if (n_sites < 0) {
    throw(std::invalid_argument("Found n_sites < 0"));
  }
}

template <typename bit_t>
combinatorics::SubsetsIterator<bit_t> BasisNoSz<bit_t>::begin() const {
  return begin_;
}

template <typename bit_t>
combinatorics::SubsetsIterator<bit_t> BasisNoSz<bit_t>::end() const {
  return end_;
}

template <typename bit_t> int64_t BasisNoSz<bit_t>::dim() const {
  return size_;
}

template <typename bit_t> int64_t BasisNoSz<bit_t>::size() const {
  return size_;
}

template <typename bit_t> int64_t BasisNoSz<bit_t>::n_sites() const {
  return n_sites_;
}

template <typename bit_t>
bool BasisNoSz<bit_t>::operator==(BasisNoSz<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_);
}

template <typename bit_t>
bool BasisNoSz<bit_t>::operator!=(BasisNoSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BasisNoSz<uint32_t>;
template class BasisNoSz<uint64_t>;

} // namespace xdiag::basis::spinhalf
