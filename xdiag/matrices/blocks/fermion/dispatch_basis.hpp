// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices {

// Concrete basis types for Fermion. Must mirror the kernel instantiation groups
// in matrices/blocks/fermion/kernels.cpp.
template <typename func_t>
void dispatch_basis(xdiag::Fermion const &block_in,
                    xdiag::Fermion const &block_out, func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;
  dispatch_basis_types<BasisOnTheFly<Subsets<uint32_t>>,           //
                       BasisOnTheFly<Subsets<uint64_t>>,           //
                       BasisOnTheFly<Combinations<uint32_t>>,      //
                       BasisOnTheFly<Combinations<uint64_t>>,      //
                       BasisOnTheFly<LinTable<uint32_t>>,          //
                       BasisOnTheFly<LinTable<uint64_t>>,          //
                       BasisOnTheFly<Combinations<BitsetDynamic>>, //
                       BasisOnTheFly<Combinations<BitsetStatic2>>, //
                       BasisOnTheFly<Combinations<BitsetStatic4>>, //
                       BasisOnTheFly<Combinations<BitsetStatic8>>, //
                       BasisSymmetric<Subsets<uint32_t>>,          //
                       BasisSymmetric<Subsets<uint64_t>>,          //
                       BasisSymmetric<Combinations<uint32_t>>,     //
                       BasisSymmetric<Combinations<uint64_t>>,     //
                       BasisSymmetric<LinTable<uint32_t>>,         //
                       BasisSymmetric<LinTable<uint64_t>>,         //
                       BasisSymmetric<Combinations<BitsetDynamic>>, //
                       BasisSymmetric<Combinations<BitsetStatic2>>, //
                       BasisSymmetric<Combinations<BitsetStatic4>>, //
                       BasisSymmetric<Combinations<BitsetStatic8>>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::matrices
