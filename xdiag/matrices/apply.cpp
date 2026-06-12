// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/matrices/blocks/boson/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/fermion/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/kernels.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// The apply function routes through two layers of type erasure before reaching
// the actual numerical kernel:
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   apply(op_t, Block, mat_t, Block, mat_t)
//     Block is a std::variant<Spinhalf, Boson, ...>. visit_same_type unwraps
//     both block_in and block_out to their concrete type (enforcing they
//     match), then calls apply_impl(OpSum, ConcreteBlock, ...). Op/Monomial are
//     first promoted to OpSum here.
//
// Layer 2 — apply_impl<block_t>: one body for every block type. The block's
//     dispatch_basis overload (matrices/blocks/<block>/dispatch_basis.hpp)
//     resolves the shared_ptr<Basis> to a concrete BasisOnTheFly<...> via a
//     type-id table. The numerical kernel is selected by block_t in
//     kernels.cpp. A block type with no dispatch_basis overload is a compile
//     error here, never a silent fallback that re-resolves to the Block
//     overload and recurses.
//
// Kernel — matrices::apply<block_t, basis_t, mat_t>(ops, basis_in, mat_in,
//          basis_out, mat_out)
//     Defined in kernels.cpp; only declared here. Each (block_t, basis_t)
//     specialisation is compiled in its own translation unit (instantiation
//     groups), keeping this file free of heavy template instantiation.

namespace xdiag {

// Layer 2 implementation — internal, called only from apply(op_t, Block, ...).
// One body for every block type: the dispatch_basis overload supplies the basis
// dispatch (the numerical kernel is selected by block_t in kernels.cpp). A block
// type with no dispatch_basis overload is a compile error here.
template <typename block_t, typename vec_t>
static void apply_impl(OpSum const &ops, block_t const &block_in,
                       vec_t const &vec_in, block_t const &block_out,
                       vec_t &vec_out) try {
  vec_out.zeros();
  // Layer 2: unwrap the basis pointer to a concrete BasisOnTheFly<...> type.
  matrices::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        // Kernel: definition is in kernels.cpp, instantiated per basis type.
        matrices::apply<block_t>(ops, basis_in, vec_in, basis_out, vec_out);
      });
}
XDIAG_CATCH

template <typename op_t, typename mat_t>
void apply(op_t const &ops, Block const &block_in, mat_t const &vec_in,
           Block const &block_out, mat_t &vec_out) try {
  // Layer 1: unwrap the Block variant (op_t is promoted to OpSum inside) and
  // forward to the block-generic apply_impl.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        apply_impl(OpSum(ops), bin, vec_in, bout, vec_out);
      },
      "Type mismatch of Block types");
}
XDIAG_CATCH

#define INSTANTIATE_XDIAG_APPLY(OP_TYPE, MAT_TYPE)                             \
  template void apply(OP_TYPE const &, Block const &, MAT_TYPE const &,        \
                      Block const &, MAT_TYPE &);

using namespace arma;
INSTANTIATE_XDIAG_APPLY(Op, vec)
INSTANTIATE_XDIAG_APPLY(Op, cx_vec)
INSTANTIATE_XDIAG_APPLY(Op, mat)
INSTANTIATE_XDIAG_APPLY(Op, cx_mat)

INSTANTIATE_XDIAG_APPLY(Monomial, vec)
INSTANTIATE_XDIAG_APPLY(Monomial, cx_vec)
INSTANTIATE_XDIAG_APPLY(Monomial, mat)
INSTANTIATE_XDIAG_APPLY(Monomial, cx_mat)

INSTANTIATE_XDIAG_APPLY(OpSum, vec)
INSTANTIATE_XDIAG_APPLY(OpSum, cx_vec)
INSTANTIATE_XDIAG_APPLY(OpSum, mat)
INSTANTIATE_XDIAG_APPLY(OpSum, cx_mat)

#undef INSTANTIATE_XDIAG_APPLY
} // namespace xdiag
