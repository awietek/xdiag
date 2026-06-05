// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebras/algebra.hpp>
#include <xdiag/matrices/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// Returns whether ops equals its Hermitian conjugate up to normal ordering with
// respect to the given algebra, i.e. isapprox(ops, hc(ops), algebra). The
// tolerance tol is forwarded to isapprox as both its relative and absolute
// tolerance.
bool ishermitian(OpSum const &ops, algebra::Algebra const &algebra,
                 double tol = 1e-12);

// Convenience overload: uses the algebra canonically associated with the block
// (see algebra(Block)).
bool ishermitian(OpSum const &ops, Block const &block, double tol = 1e-12);

// Block-aware overloads for already-assembled sparse matrices: the matrix
// stores its hermiticity flag, so the block and tolerance are ignored. These
// let templated callers that may receive either an OpSum or a sparse matrix
// invoke ishermitian(op, block) uniformly.
template <typename idx_t, typename coeff_t>
bool ishermitian(CSRMatrix<idx_t, coeff_t> const &A, Block const &block,
                 double tol = 1e-12);
template <typename idx_t, typename coeff_t>
bool ishermitian(CSCMatrix<idx_t, coeff_t> const &A, Block const &block,
                 double tol = 1e-12);
template <typename idx_t, typename coeff_t>
bool ishermitian(COOMatrix<idx_t, coeff_t> const &A, Block const &block,
                 double tol = 1e-12);

} // namespace xdiag
