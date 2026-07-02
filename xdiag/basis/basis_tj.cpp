// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_tj.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis {

// np: a single shared compressed dn fiber, linear offsets.
template <typename enumeration_t>
BasistJ<enumeration_t>::BasistJ(enumeration_t const &enum_up,
                                enumeration_t const &enum_dncs) try
    : basis_up_(enum_up) {
  dncs_.push_back(BasisOnTheFly<enumeration_t>(enum_dncs));
  int64_t n_ups = basis_up_.size();
  int64_t fiber = dncs_[0].size();
  ups_offset_.resize(n_ups);
  for (int64_t i = 0; i < n_ups; ++i) {
    ups_offset_[i] = i * fiber;
  }
  size_ = n_ups * fiber;
}
XDIAG_CATCH

// no-np: ups over all subsets; the compressed dn fiber depends only on
// popcount(ups), so store the nsites+1 distinct fibers (Subsets on the
// complement) and accumulate the prefix-sum offsets. Only valid for Subsets.
template <typename enumeration_t>
BasistJ<enumeration_t>::BasistJ(enumeration_t const &enum_up) try
    : basis_up_(enum_up) {
  if constexpr (combinatorics::is_subsets_v<enumeration_t>) {
    int64_t nsites = basis_up_.nsites();
    dncs_.reserve(nsites + 1);
    for (int64_t m = 0; m <= nsites; ++m) {
      dncs_.push_back(BasisOnTheFly<enumeration_t>(enumeration_t(nsites - m)));
    }
    ups_offset_.resize(basis_up_.size());
    int64_t off = 0;
    int64_t i = 0;
    for (bit_t ups : basis_up_) {
      ups_offset_[i] = off;
      off += dncs_[bits::popcount(ups)].size();
      ++i;
    }
    size_ = off;
  } else {
    XDIAG_THROW("single-argument BasistJ constructor is only valid for the "
                "non-number-conserving (Subsets) enumeration");
  }
}
XDIAG_CATCH

template <typename enumeration_t> int64_t BasistJ<enumeration_t>::size() const {
  return size_;
}

template <typename enumeration_t>
int64_t BasistJ<enumeration_t>::nsites() const {
  return basis_up_.nsites();
}

template <typename enumeration_t>
typename BasistJ<enumeration_t>::iterator_t
BasistJ<enumeration_t>::begin() const {
  return iterator_t(*this, true);
}

template <typename enumeration_t>
typename BasistJ<enumeration_t>::iterator_t
BasistJ<enumeration_t>::end() const {
  return iterator_t(*this, false);
}

template <typename enumeration_t>
bool BasistJ<enumeration_t>::operator==(
    BasistJ<enumeration_t> const &rhs) const {
  return (basis_up_ == rhs.basis_up_) && (dncs_ == rhs.dncs_);
}

template <typename enumeration_t>
bool BasistJ<enumeration_t>::operator!=(
    BasistJ<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename enumeration_t>
BasistJIterator<enumeration_t>::BasistJIterator(
    BasistJ<enumeration_t> const &basis, bool begin)
    : basis_(&basis), up_(basis.basis_up().begin()),
      up_end_(basis.basis_up().end()), nsites_(basis.nsites()) {
  if (begin && (basis.size() > 0)) {
    auto const &fiber = basis.basis_dncs(*up_);
    dn_ = fiber.begin();
    dn_end_ = fiber.end();
  } else {
    up_ = up_end_;
  }
}

template <typename enumeration_t>
BasistJIterator<enumeration_t> &BasistJIterator<enumeration_t>::operator++() {
  ++dn_;
  if (!(dn_ != dn_end_)) { // current dn fiber exhausted -> next ups
    ++up_;
    if (up_ != up_end_) {
      auto const &fiber = basis_->basis_dncs(*up_);
      dn_ = fiber.begin();
      dn_end_ = fiber.end();
    }
  }
  return *this;
}

template <typename enumeration_t>
std::pair<typename BasistJIterator<enumeration_t>::bit_t,
          typename BasistJIterator<enumeration_t>::bit_t>
BasistJIterator<enumeration_t>::operator*() const {
  bit_t ups = *up_;
  return {ups, tj_decompress_dns(ups, *dn_, nsites_)};
}

template <typename enumeration_t>
bool BasistJIterator<enumeration_t>::operator!=(
    BasistJIterator<enumeration_t> const &rhs) const {
  // at end (up_ == up_end_) the dn iterator is irrelevant
  return (up_ != rhs.up_) || ((up_ != up_end_) && (dn_ != rhs.dn_));
}

template <typename enumeration_t>
int64_t BasistJ<enumeration_t>::index(ProductState const &pstate) const {
  auto [ups, dns] = pair_from_pstate<bit_t>(pstate, nsites());
  bit_t dncs = tj_compress_dns(ups, dns, nsites());
  return ups_offset(index_up(ups)) + index_dncs(ups, dncs);
}

} // namespace xdiag::basis

using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_BASIS_TJ(ENUMERATION)                                      \
  template class xdiag::basis::BasistJ<ENUMERATION>;                           \
  template class xdiag::basis::BasistJIterator<ENUMERATION>;

INSTANTIATE_BASIS_TJ(LinTable<uint32_t>)
INSTANTIATE_BASIS_TJ(LinTable<uint64_t>)
INSTANTIATE_BASIS_TJ(Combinations<uint32_t>)
INSTANTIATE_BASIS_TJ(Combinations<uint64_t>)
INSTANTIATE_BASIS_TJ(Combinations<BitsetDynamic>)
INSTANTIATE_BASIS_TJ(Subsets<uint32_t>)
INSTANTIATE_BASIS_TJ(Subsets<uint64_t>)

#undef INSTANTIATE_BASIS_TJ
