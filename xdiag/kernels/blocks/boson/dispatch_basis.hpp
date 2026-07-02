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
#include <xdiag/kernels/dispatcher.hpp>

namespace xdiag::kernels {

// Concrete basis types for Boson. Must mirror the kernel instantiation groups in
// kernels/blocks/boson/kernels.cpp. Only the BitArray / BitArrayLong widths
// {1,2,3,4,8} are enumerated; widths 5-7 are promoted to 8 at dispatch time (see
// promote_nlocalbits in blocks/boson.cpp), so they are never constructed.
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
      BasisOnTheFly<SchaeferTable<BitArray8>>,
      BasisOnTheFly<BoundedPartitions<BitArray1>>,
      BasisOnTheFly<BoundedPartitions<BitArray2>>,
      BasisOnTheFly<BoundedPartitions<BitArray3>>,
      BasisOnTheFly<BoundedPartitions<BitArray4>>,
      BasisOnTheFly<BoundedPartitions<BitArray8>>,
      BasisOnTheFly<BoundedMultisets<BitArray1>>,
      BasisOnTheFly<BoundedMultisets<BitArray2>>,
      BasisOnTheFly<BoundedMultisets<BitArray3>>,
      BasisOnTheFly<BoundedMultisets<BitArray4>>,
      BasisOnTheFly<BoundedMultisets<BitArray8>>,
      BasisSymmetric<SchaeferTable<BitArray1>>,
      BasisSymmetric<SchaeferTable<BitArray2>>,
      BasisSymmetric<SchaeferTable<BitArray3>>,
      BasisSymmetric<SchaeferTable<BitArray4>>,
      BasisSymmetric<SchaeferTable<BitArray8>>,
      BasisSymmetric<BoundedPartitions<BitArray1>>,
      BasisSymmetric<BoundedPartitions<BitArray2>>,
      BasisSymmetric<BoundedPartitions<BitArray3>>,
      BasisSymmetric<BoundedPartitions<BitArray4>>,
      BasisSymmetric<BoundedPartitions<BitArray8>>,
      BasisSymmetric<BoundedMultisets<BitArray1>>,
      BasisSymmetric<BoundedMultisets<BitArray2>>,
      BasisSymmetric<BoundedMultisets<BitArray3>>,
      BasisSymmetric<BoundedMultisets<BitArray4>>,
      BasisSymmetric<BoundedMultisets<BitArray8>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong1>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong2>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong3>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong4>>,
      BasisOnTheFly<BoundedPartitions<BitArrayLong8>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong1>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong2>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong3>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong4>>,
      BasisOnTheFly<BoundedMultisets<BitArrayLong8>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong1>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong2>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong3>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong4>>,
      BasisSymmetric<BoundedPartitions<BitArrayLong8>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong1>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong2>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong3>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong4>>,
      BasisSymmetric<BoundedMultisets<BitArrayLong8>>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::kernels
