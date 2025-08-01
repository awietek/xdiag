// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <tuple>

#include <xdiag/common.hpp>

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

XDIAG_API double eigval0(OpSum const &ops, Block const &block,
                         double precision = 1e-12,
                         int64_t max_iterations = 1000,
                         int64_t random_seed = 42);

template <typename idx_t, typename coeff_t>
XDIAG_API double eigval0(CSRMatrix<idx_t, coeff_t> const &ops,
                         Block const &block, double precision = 1e-12,
                         int64_t max_iterations = 1000,
                         int64_t random_seed = 42);

std::tuple<double, State> eig0(OpSum const &ops, Block const &block,
                               double precision = 1e-12,
                               int64_t max_iterations = 1000,
                               int64_t random_seed = 42);

template <typename idx_t, typename coeff_t>
std::tuple<double, State> eig0(CSRMatrix<idx_t, coeff_t> const &ops,
                               Block const &block, double precision = 1e-12,
                               int64_t max_iterations = 1000,
                               int64_t random_seed = 42);

} // namespace xdiag
