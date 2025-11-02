// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_np.hpp"
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>

namespace xdiag::basis::electron {

using namespace combinatorics;

template <typename bit_t>
BasisNp<bit_t>::BasisNp(int nsites, int nup, int ndn) try
    : nsites_(nsites), nup_(nup), ndn_(ndn), size_ups_(binomial(nsites, nup)),
      size_dns_(binomial(nsites, ndn)), size_(size_ups_ * size_dns_),
      lintable_ups_(nsites, nup), lintable_dns_(nsites, ndn) {
  check_nsites_work_with_bits<bit_t>(nsites_);
  if (nsites < 0) {
    XDIAG_THROW("nsites < 0");
  } else if ((nup < 0) || (nup > nsites)) {
    XDIAG_THROW("Invalid value of nup");
  } else if ((ndn < 0) || (ndn > nsites)) {
    XDIAG_THROW("Invalid value of ndn");
  }
}
XDIAG_CATCH

template <typename bit_t> int BasisNp<bit_t>::nsites() const { return nsites_; }
template <typename bit_t> int BasisNp<bit_t>::nup() const { return nup_; }
template <typename bit_t> int BasisNp<bit_t>::ndn() const { return ndn_; }

template <typename bit_t> int64_t BasisNp<bit_t>::size_ups() const {
  return size_ups_;
}
template <typename bit_t> int64_t BasisNp<bit_t>::size_dns() const {
  return size_dns_;
}

template <typename bit_t> int64_t BasisNp<bit_t>::dim() const { return size_; }
template <typename bit_t> int64_t BasisNp<bit_t>::size() const { return size_; }
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::begin() const {
  return iterator_t(nsites_, nup_, ndn_, true);
}
template <typename bit_t>
typename BasisNp<bit_t>::iterator_t BasisNp<bit_t>::end() const {
  return iterator_t(nsites_, nup_, ndn_, false);
}

template <typename bit_t>
Combinations<bit_t> BasisNp<bit_t>::states_ups() const {
  return Combinations<bit_t>(nsites_, nup_);
}

template <typename bit_t>
Combinations<bit_t> BasisNp<bit_t>::states_dns() const {
  return Combinations<bit_t>(nsites_, ndn_);
}

template <typename bit_t>
CombinationsIndex<bit_t> BasisNp<bit_t>::states_indices_ups() const {
  return CombinationsIndex<bit_t>(nsites_, nup_);
}

template <typename bit_t>
CombinationsIndex<bit_t> BasisNp<bit_t>::states_indices_dns() const {
  return CombinationsIndex<bit_t>(nsites_, ndn_);
}

#ifdef _OPENMP
template <typename bit_t>
CombinationsThread<bit_t> BasisNp<bit_t>::states_ups_thread() const {
  return CombinationsThread<bit_t>(nsites_, nup_);
}

template <typename bit_t>
CombinationsThread<bit_t> BasisNp<bit_t>::states_dns_thread() const {
  return CombinationsThread<bit_t>(nsites_, ndn_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
BasisNp<bit_t>::states_indices_ups_thread() const {
  return CombinationsIndexThread<bit_t>(nsites_, nup_);
}

template <typename bit_t>
CombinationsIndexThread<bit_t>
BasisNp<bit_t>::states_indices_dns_thread() const {
  return CombinationsIndexThread<bit_t>(nsites_, ndn_);
}
#endif

template <typename bit_t>
bool BasisNp<bit_t>::operator==(BasisNp const &rhs) const {
  return (nsites_ == rhs.nsites_) && (nup_ == rhs.nup_) && (ndn_ == rhs.ndn_);
}

template <typename bit_t>
bool BasisNp<bit_t>::operator!=(BasisNp const &rhs) const {
  return !operator==(rhs);
}

template class BasisNp<uint32_t>;
template class BasisNp<uint64_t>;

template <typename bit_t>
BasisNpIterator<bit_t>::BasisNpIterator(int64_t nsites, int64_t nup,
                                        int64_t ndn, bool begin) {
  begindns_ = ((bit_t)1 << ndn) - 1;
  bit_t beginups = ((bit_t)1 << nup) - 1;
  bit_t end_ups = get_next_pattern(beginups << (nsites - nup));
  end_dns_ = get_next_pattern(begindns_ << (nsites - ndn));
  ups_ = begin ? beginups : end_ups;
  dns_ = begindns_;
}

template <typename bit_t>
BasisNpIterator<bit_t> &BasisNpIterator<bit_t>::operator++() {
  dns_ = get_next_pattern(dns_);
  if (dns_ == end_dns_) {
    dns_ = begindns_;
    ups_ = get_next_pattern(ups_);
  }
  return *this;
}

template <typename bit_t>
std::pair<bit_t, bit_t> BasisNpIterator<bit_t>::operator*() const {
  return {ups_, dns_};
}

template <typename bit_t>
bool BasisNpIterator<bit_t>::operator!=(
    BasisNpIterator<bit_t> const &rhs) const {
  return (ups_ != rhs.ups_) || (dns_ != rhs.dns_);
}

template class BasisNpIterator<uint32_t>;
template class BasisNpIterator<uint64_t>;

} // namespace xdiag::basis::electron
