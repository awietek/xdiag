// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <tuple>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Ground state calculations using Lanczos algorithm
XDIAG_API double eigval0(OpSum const &ops, Block const &block,
                         double precision = 1e-12,
                         int64_t max_iterations = 1000,
                         int64_t random_seed = 42);

template <typename idx_t, typename coeff_t>
XDIAG_API double eigval0(CSRMatrix<idx_t, coeff_t> const &ops,
                         Block const &block, double precision = 1e-12,
                         int64_t max_iterations = 1000,
                         int64_t random_seed = 42);

XDIAG_API std::tuple<double, State> eig0(OpSum const &ops, Block const &block,
                                         double precision = 1e-12,
                                         int64_t max_iterations = 1000,
                                         int64_t random_seed = 42);

template <typename idx_t, typename coeff_t>
XDIAG_API std::tuple<double, State>
eig0(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block,
     double precision = 1e-12, int64_t max_iterations = 1000,
     int64_t random_seed = 42);

// Excited state calculations using LOBPCG algorithm
XDIAG_API arma::vec eigvals(OpSum const &ops, Block const &block, int64_t neigs,
                            double precision = 1e-12,
                            int64_t max_iterations = 1000,
                            int64_t random_seed = 42);

template <typename idx_t, typename coeff_t>
XDIAG_API arma::vec
eigvals(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block, int64_t neigs,
        double precision = 1e-12, int64_t max_iterations = 1000,
        int64_t random_seed = 42);

XDIAG_API std::tuple<arma::vec, State> eigs(OpSum const &ops,
                                            Block const &block, int64_t neigs,
                                            double precision = 1e-12,
                                            int64_t max_iterations = 1000,
                                            int64_t random_seed = 42);

template <typename idx_t, typename coeff_t>
XDIAG_API std::tuple<arma::vec, State>
eigs(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block, int64_t neigs,
     double precision = 1e-12, int64_t max_iterations = 1000,
     int64_t random_seed = 42);

} // namespace xdiag
