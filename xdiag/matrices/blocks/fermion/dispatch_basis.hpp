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

namespace xdiag::matrices::fermion {

template <typename func_t>
void dispatch_basis(xdiag::Fermion const &block_in,
                    xdiag::Fermion const &block_out, func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;

  matrices::Dispatcher d;
#define ADD(B) d.add<B>([&](B const &bin, B const &bout) { fn(bin, bout); });
  ADD(BasisOnTheFly<Subsets<uint32_t>>)
  ADD(BasisOnTheFly<Subsets<uint64_t>>)
  ADD(BasisOnTheFly<Combinations<uint32_t>>)
  ADD(BasisOnTheFly<Combinations<uint64_t>>)
  ADD(BasisOnTheFly<LinTable<uint32_t>>)
  ADD(BasisOnTheFly<LinTable<uint64_t>>)
  ADD(BasisOnTheFly<Combinations<BitsetDynamic>>)
  ADD(BasisOnTheFly<Combinations<BitsetStatic2>>)
  ADD(BasisOnTheFly<Combinations<BitsetStatic4>>)
  ADD(BasisOnTheFly<Combinations<BitsetStatic8>>)
  ADD(BasisSymmetric<Subsets<uint32_t>>)
  ADD(BasisSymmetric<Subsets<uint64_t>>)
  ADD(BasisSymmetric<Combinations<uint32_t>>)
  ADD(BasisSymmetric<Combinations<uint64_t>>)
  ADD(BasisSymmetric<LinTable<uint32_t>>)
  ADD(BasisSymmetric<LinTable<uint64_t>>)
  ADD(BasisSymmetric<Combinations<BitsetDynamic>>)
  ADD(BasisSymmetric<Combinations<BitsetStatic2>>)
  ADD(BasisSymmetric<Combinations<BitsetStatic4>>)
  ADD(BasisSymmetric<Combinations<BitsetStatic8>>)
#undef ADD
  d.dispatch(block_in.basis(), block_out.basis());
}

} // namespace xdiag::matrices::fermion
