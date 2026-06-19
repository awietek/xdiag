// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_sublattice.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/blocks/spinhalf/matrix_generic.hpp>
#include <xdiag/matrices/kernels_generic.hpp>

namespace xdiag::matrices {

// The only place spinhalf::matrix_generic (and the spinhalf term headers) is
// pulled in. A change to a spinhalf term recompiles only this source's groups.
template <> struct matrix_kernel<Spinhalf> {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    spinhalf::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
  }
};

} // namespace xdiag::matrices

using namespace arma;
using namespace xdiag;
using namespace xdiag::basis;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::matrices;

//
// BasisOnTheFly Instantiations
//

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_2)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Combinations<BitsetStatic2>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_4)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Combinations<BitsetStatic4>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_8)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisOnTheFly<Combinations<BitsetStatic8>>)
// END_INSTANTIATION_GROUP

//
// BasisSymmetric Instantiations
//

// BEGIN_INSTANTIATION_GROUP(symmetric_subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_static_2)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Combinations<BitsetStatic2>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_static_4)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Combinations<BitsetStatic4>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(symmetric_combinations_bitset_static_8)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSymmetric<Combinations<BitsetStatic8>>)
// END_INSTANTIATION_GROUP

//
// BasisSublattice Instantiations
//

// BEGIN_INSTANTIATION_GROUP(sublattice_uint32_t_1)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice32<1>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint32_t_2)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice32<2>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint32_t_3)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice32<3>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint32_t_4)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice32<4>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint32_t_5)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice32<5>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint64_t_1)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice64<1>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint64_t_2)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice64<2>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint64_t_3)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice64<3>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint64_t_4)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice64<4>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(sublattice_uint64_t_5)
XDIAG_INSTANTIATE_KERNELS(Spinhalf, BasisSublattice64<5>)
// END_INSTANTIATION_GROUP
