// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "basis_onthefly.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
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
bool BasisOnTheFly<enumeration_t>::operator==(
    BasisOnTheFly<enumeration_t> const &rhs) const {
  return enumeration_ == rhs.enumeration_;
}

template <typename enumeration_t>
bool BasisOnTheFly<enumeration_t>::operator!=(
    BasisOnTheFly<enumeration_t> const &rhs) const {
  return !operator==(rhs);
}

} // namespace xdiag::basis

using namespace xdiag;
using namespace basis;
using namespace combinatorics;
using namespace bits;

#define INSTANTIATE(E) template class basis::BasisOnTheFly<E>;

// SchaeferTable / BoundedMultisets / BoundedPartitions share the BitArray
// backends, so list those once per class template.
#define INSTANTIATE_BITARRAY(Tmpl)                                             \
  INSTANTIATE(Tmpl<BitArray1>)                                                 \
  INSTANTIATE(Tmpl<BitArray2>)                                                 \
  INSTANTIATE(Tmpl<BitArray3>)                                                 \
  INSTANTIATE(Tmpl<BitArray4>)                                                 \
  INSTANTIATE(Tmpl<BitArray5>)                                                 \
  INSTANTIATE(Tmpl<BitArray6>)                                                 \
  INSTANTIATE(Tmpl<BitArray7>)                                                 \
  INSTANTIATE(Tmpl<BitArray8>)

#define INSTANTIATE_BITARRAY_LONG(Tmpl)                                        \
  INSTANTIATE(Tmpl<BitArrayLong1>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong2>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong3>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong4>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong5>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong6>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong7>)                                             \
  INSTANTIATE(Tmpl<BitArrayLong8>)

// BEGIN_INSTANTIATION_GROUP(spinhalf)
INSTANTIATE(Subsets<uint32_t>)
INSTANTIATE(Subsets<uint64_t>)
INSTANTIATE(Combinations<uint32_t>)
INSTANTIATE(Combinations<uint64_t>)
INSTANTIATE(Combinations<BitsetStatic2>)
INSTANTIATE(Combinations<BitsetStatic4>)
INSTANTIATE(Combinations<BitsetStatic8>)
INSTANTIATE(Combinations<BitsetDynamic>)
INSTANTIATE(LinTable<uint32_t>)
INSTANTIATE(LinTable<uint64_t>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(schaefer_table)
INSTANTIATE_BITARRAY(SchaeferTable)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_multisets)
INSTANTIATE_BITARRAY(BoundedMultisets)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_multisets_long)
INSTANTIATE_BITARRAY_LONG(BoundedMultisets)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_partitions)
INSTANTIATE_BITARRAY(BoundedPartitions)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_partitions_long)
INSTANTIATE_BITARRAY_LONG(BoundedPartitions)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_BITARRAY
#undef INSTANTIATE
