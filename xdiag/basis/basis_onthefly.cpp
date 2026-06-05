// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_onthefly.hpp"

#include <xdiag/basis/to_product_state.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis {
template <typename enumeration_t>
BasisOnTheFly<enumeration_t>::BasisOnTheFly(
    enumeration_t const &enumeration) try
    : enumeration_(enumeration) {}
XDIAG_CATCH

template <typename enumeration_t>
int64_t BasisOnTheFly<enumeration_t>::size() const {
  return enumeration_.size();
}

template <typename enumeration_t>
int64_t BasisOnTheFly<enumeration_t>::nsites() const {
  return enumeration_.n();
}

template <typename enumeration_t>
int64_t BasisOnTheFly<enumeration_t>::d() const {
  return enumeration_.d();
}

template <typename enumeration_t>
int64_t BasisOnTheFly<enumeration_t>::index(bit_t bits) const {
  return enumeration_.index(bits);
}

template <typename enumeration_t>
typename BasisOnTheFly<enumeration_t>::iterator_t
BasisOnTheFly<enumeration_t>::begin() const {
  return enumeration_.begin();
}

template <typename enumeration_t>
typename BasisOnTheFly<enumeration_t>::iterator_t
BasisOnTheFly<enumeration_t>::end() const {
  return enumeration_.end();
}

template <typename enumeration_t>
ProductState BasisOnTheFly<enumeration_t>::product_state(
    int64_t idx, std::vector<std::string> const &dict) const {
  return to_product_state(enumeration_, idx, dict);
}

template <typename enumeration_t>
bool BasisOnTheFly<enumeration_t>::operator==(
    BasisOnTheFly<enumeration_t> const &rhs) const {
  return enumeration_ == rhs.enumeration_;
}

template <typename enumeration_t>
bool BasisOnTheFly<enumeration_t>::operator!=(
    BasisOnTheFly<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

using namespace combinatorics;
using namespace bits;
template class BasisOnTheFly<Subsets<uint32_t>>;
template class BasisOnTheFly<Subsets<uint64_t>>;
template class BasisOnTheFly<Combinations<uint32_t>>;
template class BasisOnTheFly<Combinations<uint64_t>>;
template class BasisOnTheFly<Combinations<BitsetStatic2>>;
template class BasisOnTheFly<Combinations<BitsetStatic4>>;
template class BasisOnTheFly<Combinations<BitsetStatic8>>;
template class BasisOnTheFly<Combinations<BitsetDynamic>>;
template class BasisOnTheFly<LinTable<uint32_t>>;
template class BasisOnTheFly<LinTable<uint64_t>>;

} // namespace xdiag::basis
