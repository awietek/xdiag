// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/matrices/dispatcher.hpp>

namespace xdiag::matrices {

// Concrete basis types for Electron. Must mirror the kernel instantiation groups
// in matrices/blocks/electron/kernels.cpp and the BasisElectron instantiations
// in basis/basis_electron.cpp. Non-symmetric only (no BasisElectronSymmetric
// yet).
template <typename func_t>
void dispatch_basis(xdiag::Electron const &block_in,
                    xdiag::Electron const &block_out, func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;
  dispatch_basis_types<BasisElectron<Subsets<uint32_t>>,          //
                       BasisElectron<Subsets<uint64_t>>,          //
                       BasisElectron<Combinations<uint32_t>>,     //
                       BasisElectron<Combinations<uint64_t>>,     //
                       BasisElectron<LinTable<uint32_t>>,         //
                       BasisElectron<LinTable<uint64_t>>,         //
                       BasisElectron<Combinations<BitsetDynamic>>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::matrices
