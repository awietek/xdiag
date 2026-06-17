// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices {

// Concrete basis types for the tJ block. Must mirror the kernel instantiation
// groups in matrices/blocks/tj/kernels.cpp and the basis instantiations in
// basis/basis_tj.cpp. (No symmetric tJ basis yet.) dispatch_basis_types selects
// by the runtime basis type.
template <typename func_t>
void dispatch_basis(xdiag::tJ const &block_in, xdiag::tJ const &block_out,
                    func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;
  dispatch_basis_types<BasistJ<Subsets<uint32_t>>,           //
                       BasistJ<Subsets<uint64_t>>,           //
                       BasistJ<Combinations<uint32_t>>,      //
                       BasistJ<Combinations<uint64_t>>,      //
                       BasistJ<LinTable<uint32_t>>,          //
                       BasistJ<LinTable<uint64_t>>,          //
                       BasistJ<Combinations<BitsetDynamic>>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::matrices
