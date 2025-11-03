// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_no_sz.hpp"

#include <limits>

#include <xdiag/combinatorics/binomial.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_t>
BasisNoSz<bit_t>::BasisNoSz(int64_t nsites) try
    : nsites_(nsites), size_(pow(2, nsites)), begin_(0), end_(size_) {
  check_nsites_work_with_bits<bit_t>(nsites_);
  if (nsites < 0) {
    XDIAG_THROW("Found nsites < 0");
  }
}
XDIAG_CATCH

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

template <typename bit_t> int64_t BasisNoSz<bit_t>::nsites() const {
  return nsites_;
}

template <typename bit_t>
bool BasisNoSz<bit_t>::operator==(BasisNoSz<bit_t> const &rhs) const {
  return (nsites_ == rhs.nsites_);
}

template <typename bit_t>
bool BasisNoSz<bit_t>::operator!=(BasisNoSz<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BasisNoSz<uint32_t>;
template class BasisNoSz<uint64_t>;

} // namespace xdiag::basis::spinhalf
