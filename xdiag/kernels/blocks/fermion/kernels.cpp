// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/kernels/blocks/fermion/matrix_generic.hpp>
#include <xdiag/kernels/kernels_generic.hpp>

namespace xdiag::kernels {

// The only place fermion::matrix_generic (and the fermion term headers) is
// pulled in. A change to a fermion term recompiles only this source's groups.
template <> struct matrix_kernel<Fermion> {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    fermion::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
  }
};

} // namespace xdiag::kernels

using namespace arma;
using namespace xdiag;
using namespace xdiag::basis;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::kernels;

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_2)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Combinations<BitsetStatic2>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_4)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Combinations<BitsetStatic4>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_8)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisOnTheFly<Combinations<BitsetStatic8>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_static_2)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Combinations<BitsetStatic2>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_static_4)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Combinations<BitsetStatic4>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_static_8)
XDIAG_INSTANTIATE_KERNELS(Fermion, BasisSymmetric<Combinations<BitsetStatic8>>)
// END_INSTANTIATION_GROUP
