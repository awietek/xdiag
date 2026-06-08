// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices::boson {

template <typename func_t>
void dispatch_basis(xdiag::Boson const &block_in, xdiag::Boson const &block_out,
                    func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;

  matrices::Dispatcher d;
#define ADD(B) d.add<B>([&](B const &bin, B const &bout) { fn(bin, bout); });
  ADD(BasisOnTheFly<SchaeferTable<BitArray1>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray2>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray3>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray4>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray5>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray6>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray7>>)
  ADD(BasisOnTheFly<SchaeferTable<BitArray8>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray1>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray2>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray3>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray4>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray5>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray6>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray7>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArray8>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray1>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray2>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray3>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray4>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray5>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray6>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray7>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArray8>>)

  ADD(BasisSymmetric<SchaeferTable<BitArray1>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray2>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray3>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray4>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray5>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray6>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray7>>)
  ADD(BasisSymmetric<SchaeferTable<BitArray8>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray1>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray2>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray3>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray4>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray5>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray6>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray7>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArray8>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray1>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray2>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray3>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray4>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray5>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray6>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray7>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArray8>>)

  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong1>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong2>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong3>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong4>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong5>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong6>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong7>>)
  ADD(BasisOnTheFly<BoundedPartitions<BitArrayLong8>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong1>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong2>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong3>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong4>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong5>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong6>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong7>>)
  ADD(BasisOnTheFly<BoundedMultisets<BitArrayLong8>>)

  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong1>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong2>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong3>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong4>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong5>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong6>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong7>>)
  ADD(BasisSymmetric<BoundedPartitions<BitArrayLong8>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong1>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong2>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong3>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong4>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong5>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong6>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong7>>)
  ADD(BasisSymmetric<BoundedMultisets<BitArrayLong8>>)
#undef ADD
  d.dispatch(block_in.basis(), block_out.basis());
}

} // namespace xdiag::matrices::boson
