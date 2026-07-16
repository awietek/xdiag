// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_electron.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::basis {

template <typename enumeration_t>
BasisElectron<enumeration_t>::BasisElectron(enumeration_t const &enum_up,
                                            enumeration_t const &enum_dn) try
    : basis_up_(enum_up), basis_dn_(enum_dn) {
  if (enum_up.n() != enum_dn.n()) {
    XDIAG_THROW(fmt::format(
        "up and dn basis have different number of sites (up: {}), (dn: {})",
        enum_up.n(), enum_dn.n()));
  }
}
XDIAG_CATCH

template <typename enumeration_t>
int64_t BasisElectron<enumeration_t>::size() const {
  return basis_up_.size() * basis_dn_.size();
}

template <typename enumeration_t>
int64_t BasisElectron<enumeration_t>::nsites() const {
  return basis_up_.nsites();
}

template <typename enumeration_t>
int64_t BasisElectron<enumeration_t>::index(bit_t ups, bit_t dns) const {
  return basis_up_.index(ups) * basis_dn_.size() + basis_dn_.index(dns);
}

template <typename enumeration_t>
int64_t
BasisElectron<enumeration_t>::index(ProductState const &pstate) const {
  auto [ups, dns] = pair_from_pstate<bit_t>(pstate, nsites());
  return index(ups, dns);
}

template <typename enumeration_t>
typename BasisElectron<enumeration_t>::iterator_t
BasisElectron<enumeration_t>::begin() const {
  return iterator_t(basis_up_.begin(), basis_up_.end(), basis_dn_.begin(),
                    basis_dn_.begin(), basis_dn_.end());
}

template <typename enumeration_t>
typename BasisElectron<enumeration_t>::iterator_t
BasisElectron<enumeration_t>::end() const {
  return iterator_t(basis_up_.end(), basis_up_.end(), basis_dn_.begin(),
                    basis_dn_.begin(), basis_dn_.end());
}

template <typename enumeration_t>
BasisOnTheFly<enumeration_t> const &
BasisElectron<enumeration_t>::basis_up() const {
  return basis_up_;
}

template <typename enumeration_t>
BasisOnTheFly<enumeration_t> const &
BasisElectron<enumeration_t>::basis_dn() const {
  return basis_dn_;
}

template <typename enumeration_t>
bool BasisElectron<enumeration_t>::operator==(
    BasisElectron<enumeration_t> const &rhs) const {
  return (basis_up_ == rhs.basis_up_) && (basis_dn_ == rhs.basis_dn_);
}

template <typename enumeration_t>
bool BasisElectron<enumeration_t>::operator!=(
    BasisElectron<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

} // namespace xdiag::basis

using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_BASIS_ELECTRON(ENUMERATION)                                \
  template class xdiag::basis::BasisElectron<ENUMERATION>;

INSTANTIATE_BASIS_ELECTRON(Subsets<uint32_t>)
INSTANTIATE_BASIS_ELECTRON(Subsets<uint64_t>)
INSTANTIATE_BASIS_ELECTRON(LinTable<uint32_t>)
INSTANTIATE_BASIS_ELECTRON(LinTable<uint64_t>)
INSTANTIATE_BASIS_ELECTRON(Combinations<uint32_t>)
INSTANTIATE_BASIS_ELECTRON(Combinations<uint64_t>)
INSTANTIATE_BASIS_ELECTRON(Combinations<BitsetDynamic>)

#undef INSTANTIATE_BASIS_ELECTRON
