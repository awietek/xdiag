// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/kernels/dispatcher.hpp>

namespace xdiag::kernels {

// Concrete basis types for Electron. Must mirror the kernel instantiation groups
// in kernels/blocks/electron/kernels.cpp and the basis instantiations in
// basis/basis_electron.cpp and basis/basis_electron_symmetric.cpp. The
// symmetric and non-symmetric types coexist; dispatch_basis_types selects by the
// runtime basis type.
template <typename func_t>
void dispatch_basis(xdiag::Electron const &block_in,
                    xdiag::Electron const &block_out, func_t fn) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;
  dispatch_basis_types<BasisElectron<Subsets<uint32_t>>,                   //
                       BasisElectron<Subsets<uint64_t>>,                   //
                       BasisElectron<Combinations<uint32_t>>,              //
                       BasisElectron<Combinations<uint64_t>>,              //
                       BasisElectron<LinTable<uint32_t>>,                  //
                       BasisElectron<LinTable<uint64_t>>,                  //
                       BasisElectron<Combinations<BitsetDynamic>>,         //
                       BasisElectronSymmetric<Subsets<uint32_t>>,          //
                       BasisElectronSymmetric<Subsets<uint64_t>>,          //
                       BasisElectronSymmetric<Combinations<uint32_t>>,     //
                       BasisElectronSymmetric<Combinations<uint64_t>>,     //
                       BasisElectronSymmetric<LinTable<uint32_t>>,         //
                       BasisElectronSymmetric<LinTable<uint64_t>>,         //
                       BasisElectronSymmetric<Combinations<BitsetDynamic>>>(
      *block_in.basis(), *block_out.basis(), fn);
}

} // namespace xdiag::kernels
