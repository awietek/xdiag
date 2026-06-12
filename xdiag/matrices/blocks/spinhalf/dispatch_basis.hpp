// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_sublattice.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices {

// Concrete basis types for Spinhalf. Must mirror the kernel instantiation groups
// in matrices/blocks/spinhalf/kernels.cpp.
template <typename func_t>
void dispatch_basis(xdiag::Spinhalf const &block_in,
                    xdiag::Spinhalf const &block_out, func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;
  dispatch_basis_types<BasisOnTheFly<Subsets<uint32_t>>,            //
                       BasisOnTheFly<Subsets<uint64_t>>,            //
                       BasisOnTheFly<Combinations<uint32_t>>,       //
                       BasisOnTheFly<Combinations<uint64_t>>,       //
                       BasisOnTheFly<LinTable<uint32_t>>,           //
                       BasisOnTheFly<LinTable<uint64_t>>,           //
                       BasisOnTheFly<Combinations<BitsetDynamic>>,  //
                       BasisOnTheFly<Combinations<BitsetStatic2>>,  //
                       BasisOnTheFly<Combinations<BitsetStatic4>>,  //
                       BasisOnTheFly<Combinations<BitsetStatic8>>,  //
                       BasisSymmetric<Subsets<uint32_t>>,           //
                       BasisSymmetric<Subsets<uint64_t>>,           //
                       BasisSymmetric<Combinations<uint32_t>>,      //
                       BasisSymmetric<Combinations<uint64_t>>,      //
                       BasisSymmetric<LinTable<uint32_t>>,          //
                       BasisSymmetric<LinTable<uint64_t>>,          //
                       BasisSymmetric<Combinations<BitsetDynamic>>, //
                       BasisSymmetric<Combinations<BitsetStatic2>>, //
                       BasisSymmetric<Combinations<BitsetStatic4>>, //
                       BasisSymmetric<Combinations<BitsetStatic8>>, //
                       BasisSublattice32<1>, BasisSublattice32<2>,  //
                       BasisSublattice32<3>, BasisSublattice32<4>,  //
                       BasisSublattice32<5>, BasisSublattice64<1>,  //
                       BasisSublattice64<2>, BasisSublattice64<3>,  //
                       BasisSublattice64<4>, BasisSublattice64<5>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::matrices
