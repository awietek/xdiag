// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/blocks/tj/matrix_generic.hpp>
#include <xdiag/matrices/kernels_generic.hpp>

namespace xdiag::matrices {

// The only place tj::matrix_generic (and the tJ term headers) is pulled in. A
// change to a tJ term recompiles only this source's groups.
template <> struct matrix_kernel<tJ> {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    tj::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
  }
};

} // namespace xdiag::matrices

using namespace arma;
using namespace xdiag;
using namespace xdiag::basis;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::matrices;

// BEGIN_INSTANTIATION_GROUP(tj_subsets_uint32_t)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(tj_subsets_uint64_t)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(tj_combinations_uint32_t)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(tj_combinations_uint64_t)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(tj_lintable_uint32_t)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(tj_lintable_uint64_t)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(tj_combinations_bitset_dynamic)
XDIAG_INSTANTIATE_KERNELS(tJ, BasistJ<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP
