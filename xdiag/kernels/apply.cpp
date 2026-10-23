// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/blocks/dispatch_bases.hpp>
#include <xdiag/kernels/blocks/distributed/apply_distributed.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// apply(ops, block_in, vec_in, block_out, vec_out) compiles ops
// (kernels::compile: sector check, transpose, normal ordering, see compile.hpp) and
// calls apply_compiled. Iterative algorithms compile once themselves and call
// kernels::apply(CompiledOpSum, ...) on every iteration.
//
// apply_compiled unwraps the Block variant (visit_same_type), apply_block
// resolves the concrete basis types (dispatch_basis) and calls the kernel
// kernels::apply<block_t> (kernels.hpp, instantiated per basis in
// kernels/blocks/<block>/kernels.cpp). The kernels enumerate the rows
// (block_out) and write each output element from one thread (gather fill), so
// the result is bitwise reproducible for any number of threads.

namespace xdiag::kernels {

// The basis-level kernel: defined in kernels.hpp and explicitly instantiated
// per basis type in kernels/blocks/<block>/kernels.cpp, so this file never
// compiles the term kernels itself.
template <typename block_t, typename basis_t, typename mat_t>
void apply(OpSum const &ops, basis_t const &rows, mat_t const &vec_in,
           basis_t const &cols, mat_t &vec_out);

} // namespace xdiag::kernels

namespace xdiag {

// Applies the compiled ops on a concrete block type: resolves the basis types
// and calls the kernel (or the MPI path for a distributed block).
template <typename block_t, typename vec_t>
static void apply_block(OpSum const &ops, block_t const &block_in,
                        vec_t const &vec_in, block_t const &block_out,
                        vec_t &vec_out) try {
  vec_out.zeros();

  if constexpr (is_distributed_v<block_t>) {
#ifdef XDIAG_DISTRIBUTED
    // Distributed blocks use a separate matrix-free, MPI-aware apply path
    // (kernels/blocks/distributed/<block>/kernels.cpp) that operates on a
    // single column. A multi-column block (e.g. a LOBPCG search space) is
    // applied one column at a time. For distributed blocks ops is not compiled
    // (see compile.hpp): their kernels normal-order it and apply it as is.
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
    // ops is the compiled transpose(ops) (see compile.hpp): the kernels run it
    // from
    // the rows (block_out) to the columns (block_in), and the gather fill turns
    // that into vec_out = ops * vec_in.
    kernels::dispatch_basis(
        block_out, block_in, [&](auto const &rows, auto const &cols) {
          kernels::apply<block_t>(ops, rows, vec_in, cols, vec_out);
        });
  }
}
XDIAG_CATCH

// Unwraps the Block variants of the compiled ops and applies it.
template <typename mat_t>
static void apply_compiled(kernels::CompiledOpSum const &ops,
                           mat_t const &vec_in, mat_t &vec_out) try {
  utils::visit_same_type(
      ops.block_in, ops.block_out,
      [&](auto const &bin, auto const &bout) {
        apply_block(ops.ops, bin, vec_in, bout, vec_out);
      },
      "Type mismatch of Block types");
}
XDIAG_CATCH

// The public overloads, one per operator type (Op, Monomial, OpSum) and vector
// / matrix type: non-templated for the Julia wrapper generator.
#define XDIAG_DEFINE_APPLY(OP_TYPE, MAT_TYPE)                                  \
  void apply(OP_TYPE const &ops, Block const &block_in,                        \
             MAT_TYPE const &vec_in, Block const &block_out,                   \
             MAT_TYPE &vec_out) try {                                          \
    apply_compiled(kernels::compile(OpSum(ops), block_in, block_out), vec_in,  \
                   vec_out);                                                   \
  }                                                                            \
  XDIAG_CATCH

#define XDIAG_DEFINE_APPLY_ALL(OP_TYPE)                                        \
  XDIAG_DEFINE_APPLY(OP_TYPE, arma::vec)                                       \
  XDIAG_DEFINE_APPLY(OP_TYPE, arma::cx_vec)                                    \
  XDIAG_DEFINE_APPLY(OP_TYPE, arma::mat)                                       \
  XDIAG_DEFINE_APPLY(OP_TYPE, arma::cx_mat)

XDIAG_DEFINE_APPLY_ALL(Op)
XDIAG_DEFINE_APPLY_ALL(Monomial)
XDIAG_DEFINE_APPLY_ALL(OpSum)

#undef XDIAG_DEFINE_APPLY_ALL
#undef XDIAG_DEFINE_APPLY

namespace kernels {

void apply(CompiledOpSum const &ops, arma::vec const &vec_in,
           arma::vec &vec_out) try {
  apply_compiled(ops, vec_in, vec_out);
}
XDIAG_CATCH
void apply(CompiledOpSum const &ops, arma::cx_vec const &vec_in,
           arma::cx_vec &vec_out) try {
  apply_compiled(ops, vec_in, vec_out);
}
XDIAG_CATCH
void apply(CompiledOpSum const &ops, arma::mat const &vec_in,
           arma::mat &vec_out) try {
  apply_compiled(ops, vec_in, vec_out);
}
XDIAG_CATCH
void apply(CompiledOpSum const &ops, arma::cx_mat const &vec_in,
           arma::cx_mat &vec_out) try {
  apply_compiled(ops, vec_in, vec_out);
}
XDIAG_CATCH

} // namespace kernels
} // namespace xdiag
