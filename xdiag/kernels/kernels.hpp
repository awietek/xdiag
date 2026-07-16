// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/operators/opsum.hpp>

// Generic kernel declarations.  Each function is templated on the block type
// (Spinhalf, Boson, Fermion, ...). The block selects the numerical kernel:
// kernels.cpp maps block_t to that block's matrix_generic implementation via
// the internal matrix_kernel<block_t> trait.
//
// These are only DECLARED here; consumers (apply.cpp, matrix.cpp, the sparse
// builders) name them with an incomplete block type and forward to the extern
// instantiation, so they never pull in matrix_generic or the term headers.
// Definitions live in kernels.cpp; explicit instantiations are provided there
// for each (block_t, basis_t) pair using the instantiation-group mechanism.

namespace xdiag::kernels {

template <typename block_t, typename basis_t, typename mat_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out);

template <typename block_t, typename coeff_t, typename basis_t>
void matrix(OpSum const &ops, basis_t const &basis_in, basis_t const &basis_out,
            coeff_t *mat);

#ifdef _OPENMP
template <typename block_t, typename coeff_t, typename basis_t>
std::vector<int64_t> coo_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                                    basis_t const &basis_out);
#else
template <typename block_t, typename coeff_t, typename basis_t>
int64_t coo_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                       basis_t const &basis_out);
#endif

#ifdef _OPENMP
template <typename block_t, typename coeff_t, typename basis_t,
          typename idx_t>
void coo_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out,
                     std::vector<int64_t> const &nnz_thread, idx_t *rows,
                     idx_t *cols, coeff_t *data, idx_t i0);
#else
template <typename block_t, typename coeff_t, typename basis_t,
          typename idx_t>
void coo_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out, idx_t *rows, idx_t *cols,
                     coeff_t *data, idx_t i0);
#endif

template <typename block_t, typename coeff_t, typename basis_t>
std::vector<int64_t> csr_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                                    basis_t const &basis_out,
                                    bool transpose = false);

template <typename block_t, typename coeff_t, typename basis_t,
          typename idx_t>
void csr_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out, std::vector<int64_t> &offset,
                     idx_t *col, coeff_t *data, idx_t i0,
                     bool transpose = false);

} // namespace xdiag::kernels
