// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/blocks/boson/matrix_generic.hpp>
#include <xdiag/matrices/kernels_generic.hpp>

namespace xdiag::matrices {

// The only place boson::matrix_generic (and the boson term headers) is pulled
// in. A change to a boson term recompiles only this source's groups.
template <> struct matrix_kernel<Boson> {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    boson::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
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
//
// Boson Instantiations (must mirror matrices/blocks/boson/dispatch_basis.hpp)
//
//

// BEGIN_INSTANTIATION_GROUP(boson_onthefly_schaefer_table)
INSTANTIATE_KERNELS_BITARRAY(Boson, BasisOnTheFly, SchaeferTable)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_onthefly_bounded_partitions)
INSTANTIATE_KERNELS_BITARRAY(Boson, BasisOnTheFly, BoundedPartitions)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_onthefly_bounded_multisets)
INSTANTIATE_KERNELS_BITARRAY(Boson, BasisOnTheFly, BoundedMultisets)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_symmetric_schaefer_table)
INSTANTIATE_KERNELS_BITARRAY(Boson, BasisSymmetric, SchaeferTable)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_symmetric_bounded_partitions)
INSTANTIATE_KERNELS_BITARRAY(Boson, BasisSymmetric, BoundedPartitions)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_symmetric_bounded_multisets)
INSTANTIATE_KERNELS_BITARRAY(Boson, BasisSymmetric, BoundedMultisets)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_onthefly_bounded_partitions_long)
INSTANTIATE_KERNELS_BITARRAY_LONG(Boson, BasisOnTheFly, BoundedPartitions)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_onthefly_bounded_multisets_long)
INSTANTIATE_KERNELS_BITARRAY_LONG(Boson, BasisOnTheFly, BoundedMultisets)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_symmetric_bounded_partitions_long)
INSTANTIATE_KERNELS_BITARRAY_LONG(Boson, BasisSymmetric, BoundedPartitions)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(boson_symmetric_bounded_multisets_long)
INSTANTIATE_KERNELS_BITARRAY_LONG(Boson, BasisSymmetric, BoundedMultisets)
// END_INSTANTIATION_GROUP
