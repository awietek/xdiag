// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Computes vec_out = ops * vec_in where vec_in/vec_out are arma::Col or
// arma::Mat (real or complex).
// vec_out is zeroed before accumulation. block_in and block_out must be the
// same Block type.

// Comment: this is not templated over the vec/mat type, to make it
// automatically accessible to the Julia wrapper generator
XDIAG_API void apply(Op const &op, Block const &block_in,
                     arma::vec const &vec_in, Block const &block_out,
                     arma::vec &vec_out);
XDIAG_API void apply(Op const &op, Block const &block_in,
                     arma::cx_vec const &vec_in, Block const &block_out,
                     arma::cx_vec &vec_out);
XDIAG_API void apply(Op const &op, Block const &block_in,
                     arma::mat const &vec_in, Block const &block_out,
                     arma::mat &vec_out);
XDIAG_API void apply(Op const &op, Block const &block_in,
                     arma::cx_mat const &vec_in, Block const &block_out,
                     arma::cx_mat &vec_out);

XDIAG_API void apply(Monomial const &op, Block const &block_in,
                     arma::vec const &vec_in, Block const &block_out,
                     arma::vec &vec_out);
XDIAG_API void apply(Monomial const &op, Block const &block_in,
                     arma::cx_vec const &vec_in, Block const &block_out,
                     arma::cx_vec &vec_out);
XDIAG_API void apply(Monomial const &op, Block const &block_in,
                     arma::mat const &vec_in, Block const &block_out,
                     arma::mat &vec_out);
XDIAG_API void apply(Monomial const &op, Block const &block_in,
                     arma::cx_mat const &vec_in, Block const &block_out,
                     arma::cx_mat &vec_out);

XDIAG_API void apply(OpSum const &op, Block const &block_in,
                     arma::vec const &vec_in, Block const &block_out,
                     arma::vec &vec_out);
XDIAG_API void apply(OpSum const &op, Block const &block_in,
                     arma::cx_vec const &vec_in, Block const &block_out,
                     arma::cx_vec &vec_out);
XDIAG_API void apply(OpSum const &op, Block const &block_in,
                     arma::mat const &vec_in, Block const &block_out,
                     arma::mat &vec_out);
XDIAG_API void apply(OpSum const &op, Block const &block_in,
                     arma::cx_mat const &vec_in, Block const &block_out,
                     arma::cx_mat &vec_out);

namespace kernels {

// Internal (not part of the public API).
//
// Computes vec_out = ops * vec_in for a HERMITIAN ops, bitwise reproducible
// for a fixed number of OpenMP threads: the gather fill writes every output
// element from exactly one thread in a fixed order, without atomics
// (see fill_apply<true> in fill_functions.hpp).
//
// Precondition: ops must be Hermitian on the block. The caller is responsible
// for verifying this, e.g. with ishermitian(ops, block). Internally the gather
// fill computes hc(ops) * vec_in, so for a non-Hermitian ops the result is the
// adjoint applied, not ops itself, and no error is raised.
//
// Intended for iterative algorithms which regenerate a Krylov sequence and
// therefore need identical matrix-vector products across passes (Lanczos,
// LOBPCG, time evolution).

void apply_hermitian(OpSum const &ops, Block const &block_in,
                     arma::vec const &vec_in, Block const &block_out,
                     arma::vec &vec_out);
void apply_hermitian(OpSum const &ops, Block const &block_in,
                     arma::cx_vec const &vec_in, Block const &block_out,
                     arma::cx_vec &vec_out);
void apply_hermitian(OpSum const &ops, Block const &block_in,
                     arma::mat const &vec_in, Block const &block_out,
                     arma::mat &vec_out);
void apply_hermitian(OpSum const &ops, Block const &block_in,
                     arma::cx_mat const &vec_in, Block const &block_out,
                     arma::cx_mat &vec_out);

} // namespace kernels
} // namespace xdiag
