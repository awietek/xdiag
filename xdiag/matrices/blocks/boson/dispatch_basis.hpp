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

namespace xdiag::matrices {

// Concrete basis types for Boson. Must mirror the kernel instantiation groups in
// matrices/blocks/boson/kernels.cpp. Both the BitArray<1..8> short backends and
// the BitArrayLong<1..8> backends are enumerated.
template <typename func_t>
void dispatch_basis(xdiag::Boson const &block_in, xdiag::Boson const &block_out,
                    func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;
  dispatch_basis_types<
      BasisOnTheFly<SchaeferTable<BitArray1>>,
      BasisOnTheFly<SchaeferTable<BitArray2>>,
      BasisOnTheFly<SchaeferTable<BitArray3>>,
      BasisOnTheFly<SchaeferTable<BitArray4>>,
      BasisOnTheFly<SchaeferTable<BitArray5>>,
      BasisOnTheFly<SchaeferTable<BitArray6>>,
      BasisOnTheFly<SchaeferTable<BitArray7>>,
      BasisOnTheFly<SchaeferTable<BitArray8>>,
      BasisOnTheFly<BoundedPartitions<BitArray1>>,
      BasisOnTheFly<BoundedPartitions<BitArray2>>,
      BasisOnTheFly<BoundedPartitions<BitArray3>>,
      BasisOnTheFly<BoundedPartitions<BitArray4>>,
      BasisOnTheFly<BoundedPartitions<BitArray5>>,
      BasisOnTheFly<BoundedPartitions<BitArray6>>,
      BasisOnTheFly<BoundedPartitions<BitArray7>>,
      BasisOnTheFly<BoundedPartitions<BitArray8>>,
      BasisOnTheFly<BoundedMultisets<BitArray1>>,
      BasisOnTheFly<BoundedMultisets<BitArray2>>,
      BasisOnTheFly<BoundedMultisets<BitArray3>>,
      BasisOnTheFly<BoundedMultisets<BitArray4>>,
      BasisOnTheFly<BoundedMultisets<BitArray5>>,
      BasisOnTheFly<BoundedMultisets<BitArray6>>,
      BasisOnTheFly<BoundedMultisets<BitArray7>>,
      BasisOnTheFly<BoundedMultisets<BitArray8>>,
      BasisSymmetric<SchaeferTable<BitArray1>>,
      BasisSymmetric<SchaeferTable<BitArray2>>,
      BasisSymmetric<SchaeferTable<BitArray3>>,
      BasisSymmetric<SchaeferTable<BitArray4>>,
      BasisSymmetric<SchaeferTable<BitArray5>>,
      BasisSymmetric<SchaeferTable<BitArray6>>,
      BasisSymmetric<SchaeferTable<BitArray7>>,
      BasisSymmetric<SchaeferTable<BitArray8>>,
      BasisSymmetric<BoundedPartitions<BitArray1>>,
      BasisSymmetric<BoundedPartitions<BitArray2>>,
      BasisSymmetric<BoundedPartitions<BitArray3>>,
      BasisSymmetric<BoundedPartitions<BitArray4>>,
      BasisSymmetric<BoundedPartitions<BitArray5>>,
      BasisSymmetric<BoundedPartitions<BitArray6>>,
      BasisSymmetric<BoundedPartitions<BitArray7>>,
      BasisSymmetric<BoundedPartitions<BitArray8>>,
      BasisSymmetric<BoundedMultisets<BitArray1>>,
      BasisSymmetric<BoundedMultisets<BitArray2>>,
      BasisSymmetric<BoundedMultisets<BitArray3>>,
      BasisSymmetric<BoundedMultisets<BitArray4>>,
      BasisSymmetric<BoundedMultisets<BitArray5>>,
      BasisSymmetric<BoundedMultisets<BitArray6>>,
      BasisSymmetric<BoundedMultisets<BitArray7>>,
      BasisSymmetric<BoundedMultisets<BitArray8>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong1>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong2>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong3>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong4>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong5>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong6>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong7>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong8>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong1>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong2>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong3>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong4>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong5>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong6>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong7>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong8>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong1>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong2>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong3>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong4>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong5>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong6>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong7>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong8>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong1>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong2>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong3>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong4>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong5>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong6>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong7>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong8>>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::matrices
