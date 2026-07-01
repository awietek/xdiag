// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/blocks/dispatch_bases.hpp>
#include <xdiag/kernels/blocks/distributed/apply_distributed.hpp>
#include <xdiag/kernels/kernels.hpp>
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
//     dispatch_basis overload (kernels/blocks/<block>/dispatch_basis.hpp)
//     resolves the shared_ptr<Basis> to a concrete BasisOnTheFly<...> via a
//     type-id table. The numerical kernel is selected by block_t in
//     kernels.cpp. A block type with no dispatch_basis overload is a compile
//     error here, never a silent fallback that re-resolves to the Block
//     overload and recurses.
//
// Kernel — kernels::apply<block_t, basis_t, mat_t>(ops, basis_in, mat_in,
//          basis_out, mat_out)
//     Defined in kernels.cpp; only declared here. Each (block_t, basis_t)
//     specialisation is compiled in its own translation unit (instantiation
//     groups), keeping this file free of heavy template instantiation.

namespace xdiag {

// Layer 2 implementation — internal, called only from apply(op_t, Block, ...).
// One body for every block type: the dispatch_basis overload supplies the basis
// dispatch (the numerical kernel is selected by block_t in kernels.cpp). A
// block type with no dispatch_basis overload is a compile error here.
template <typename block_t, typename vec_t>
static void apply_impl(OpSum const &ops, block_t const &block_in,
                       vec_t const &vec_in, block_t const &block_out,
                       vec_t &vec_out) try {
  vec_out.zeros();

  if constexpr (is_distributed_v<block_t>) {
#ifdef XDIAG_DISTRIBUTED
    // Distributed blocks use a separate matrix-free, MPI-aware apply path
    // (kernels/blocks/distributed/<block>/kernels.cpp) that operates on a
    // single column. A multi-column block (e.g. a LOBPCG search space) is
    // applied one column at a time.
    using coeff_t = typename vec_t::elem_type;
    if constexpr (std::is_same_v<vec_t, arma::Col<coeff_t>>) {
      kernels::apply_distributed(ops, block_in, vec_in, block_out, vec_out);
    } else {
      for (arma::uword c = 0; c < vec_in.n_cols; ++c) {
        arma::Col<coeff_t> in_col = vec_in.col(c);
        arma::Col<coeff_t> out_col(vec_out.n_rows, arma::fill::zeros);
        kernels::apply_distributed(ops, block_in, in_col, block_out, out_col);
        vec_out.col(c) = out_col;
      }
    }
#endif
  } else {
    // Layer 2: unwrap the basis pointer to a concrete BasisOnTheFly<...> type.
    kernels::dispatch_basis(
        block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
          // Kernel: definition is in kernels.cpp, instantiated per basis type.
          kernels::apply<block_t>(ops, basis_in, vec_in, basis_out, vec_out);
        });
  }
}
XDIAG_CATCH

template <typename mat_t>
XDIAG_API void apply(Op const &op, Block const &block_in, mat_t const &vec_in,
                     Block const &block_out, mat_t &vec_out) {
  apply(OpSum(op), block_in, vec_in, block_out, vec_out);
}

template <typename mat_t>
XDIAG_API void apply(Monomial const &mono, Block const &block_in,
                     mat_t const &vec_in, Block const &block_out,
                     mat_t &vec_out) {
  apply(OpSum(mono), block_in, vec_in, block_out, vec_out);
}

template <typename mat_t>
void apply(OpSum const &ops, Block const &block_in, mat_t const &vec_in,
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
