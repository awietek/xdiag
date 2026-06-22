// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/kernels/blocks/electron/matrix_generic.hpp>
#include <xdiag/kernels/kernels_generic.hpp>

namespace xdiag::kernels {

// The only place electron::matrix_generic (and the electron term headers) is
// pulled in. A change to an electron term recompiles only this source's groups.
template <> struct matrix_kernel<Electron> {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    electron::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
  }
};

} // namespace xdiag::kernels

using namespace arma;
using namespace xdiag;
using namespace xdiag::basis;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::kernels;

// BEGIN_INSTANTIATION_GROUP(subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectron<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectronSymmetric<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectronSymmetric<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Electron,
                          BasisElectronSymmetric<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Electron,
                          BasisElectronSymmetric<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectronSymmetric<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Electron, BasisElectronSymmetric<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(Electron,
                          BasisElectronSymmetric<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP
