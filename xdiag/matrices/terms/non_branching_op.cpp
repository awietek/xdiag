// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "non_branching_op.hpp"

#include <tuple>
#include <vector>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::matrices {

template <typename bit_t, typename coeff_t>
NonBranchingOp<bit_t, coeff_t>::NonBranchingOp(
    std::vector<int64_t> const &sites, arma::Mat<coeff_t> const &matrix,
    int64_t d, double precision)
    : sites_(sites), d_(d), diagonal_(true) {

  int64_t dim = math::ipow(d, (int64_t)sites.size());
  if ((int64_t)matrix.n_rows != dim || (int64_t)matrix.n_cols != dim) {
    XDIAG_THROW(fmt::format(
        "Invalid matrix dimension for non-branching Op. Expected dim={}, "
        "got n_rows={}, n_cols={}.",
        dim, matrix.n_rows, matrix.n_cols));
  }

  // Build hops_, check non-branching, detect diagonality in a single pass.
  hops_.assign(dim, {int64_t(0), coeff_t(0)});
  std::vector<bool> row_used(dim, false);
  for (int64_t in = 0; in < dim; ++in) {
    bool found = false;
    for (int64_t out = 0; out < dim; ++out) {
      if (std::abs(matrix(out, in)) > precision) {
        if (found) {
          XDIAG_THROW(fmt::format(
              "Matrix is branching: column {} has multiple non-zero entries",
              in));
        }
        if (row_used[out]) {
          XDIAG_THROW(fmt::format(
              "Matrix is branching: row {} has multiple non-zero entries",
              out));
        }
        found = true;
        row_used[out] = true;
        hops_[in] = {out, matrix(out, in)};
        if (out != in) {
          diagonal_ = false;
        }
      }
    }
  }
}

template <typename bit_t, typename coeff_t>
int64_t NonBranchingOp<bit_t, coeff_t>::extract(bit_t state) const {
  int64_t local = 0;
  for (int64_t i = (int64_t)sites_.size() - 1; i >= 0; --i) {
    local = local * d_ + bits::get(state, sites_[i]);
  }
  return local;
}

template <typename bit_t, typename coeff_t>
bit_t NonBranchingOp<bit_t, coeff_t>::deposit(int64_t local,
                                              bit_t state) const {
  for (int64_t i = 0; i < (int64_t)sites_.size(); ++i) {
    bits::set(state, sites_[i], local % d_);
    local /= d_;
  }
  return state;
}

template <typename coeff_t>
static std::vector<arma::Mat<coeff_t>>
decompose_matrix_to_nonbranching(arma::Mat<coeff_t> const &mat,
                                 double precision) try {
  int64_t m = (int64_t)mat.n_rows;
  int64_t n = (int64_t)mat.n_cols;
  if (m != n) {
    XDIAG_THROW("Error: Op matrix is not square");
  }

  std::vector<std::tuple<int64_t, int64_t, coeff_t>> all_entries;

  // Get diagonal elements first so the first non-branching matrix tends to
  // contain the entire diagonal (so the diagonal fast-path in term_matrix
  // triggers cleanly).
  for (int64_t i = 0; i < n; ++i) {
    if (std::abs(mat(i, i)) > precision) {
      all_entries.push_back({i, i, mat(i, i)});
    }
  }

  // Get offdiagonal elements
  for (int64_t n_diag = 1; n_diag < n; ++n_diag) {
    for (int64_t i = 0; i < n - n_diag; ++i) {
      if (std::abs(mat(i, i + n_diag)) > precision) {
        all_entries.push_back({i, i + n_diag, mat(i, i + n_diag)});
      }
      if (std::abs(mat(i + n_diag, i)) > precision) {
        all_entries.push_back({i + n_diag, i, mat(i + n_diag, i)});
      }
    }
  }

  // Reduce to minimal number of non-branching terms via greedy bipartite
  // matching: each pass picks the maximum independent set of entries (no
  // shared row, no shared column), produces one non-branching matrix.
  std::vector<arma::Mat<coeff_t>> mats_nb;
  std::vector<bool> consumed(all_entries.size(), false);
  std::vector<bool> row_taken(n);
  std::vector<bool> col_taken(n);
  int64_t remaining = (int64_t)all_entries.size();

  while (remaining > 0) {
    std::fill(row_taken.begin(), row_taken.end(), false);
    std::fill(col_taken.begin(), col_taken.end(), false);
    arma::Mat<coeff_t> mat_nb(m, n, arma::fill::zeros);

    for (int64_t idx = 0; idx < (int64_t)all_entries.size(); ++idx) {
      if (consumed[idx]) {
        continue;
      }
      auto [row, col, coeff] = all_entries[idx];
      if (!row_taken[row] && !col_taken[col]) {
        row_taken[row] = true;
        col_taken[col] = true;
        mat_nb(row, col) = coeff;
        consumed[idx] = true;
        --remaining;
      }
    }
    mats_nb.push_back(mat_nb);
  }
  return mats_nb;
} XDIAG_CATCH

static std::vector<Matrix>
decompose_matrix_to_nonbranching(Matrix const &mat, double precision) try {
  std::vector<Matrix> mats;
  if (mat.isreal()) {
    std::vector<arma::mat> mats2 =
        decompose_matrix_to_nonbranching(mat.as<arma::mat>(), precision);
    for (auto const &mat : mats2) {
      mats.push_back(Matrix(mat));
    }
  } else {
    std::vector<arma::cx_mat> mats2 =
        decompose_matrix_to_nonbranching(mat.as<arma::cx_mat>(), precision);
    for (auto const &mat : mats2) {
      mats.push_back(Matrix(mat));
    }
  }
  return mats;
} XDIAG_CATCH

template <typename bit_t, typename coeff_t>
std::vector<NonBranchingOp<bit_t, coeff_t>>
non_branching_ops(Coeff const &cpl, Op const &op, int64_t d,
                  double precision) try {
  if (cpl.isstring()) {
    XDIAG_THROW("Cannot convert Op to NonBranchingOps. Coeff found to be a "
                "string. To convert it to a NonBranchingOp the coupling must "
                "be either a real or complex number.");
  }

  std::vector<NonBranchingOp<bit_t, coeff_t>> ops;
  if (op.hasmatrix()) {
    auto mat = op.matrix() * cpl.scalar();
    auto sites = op.sites();
    auto mats_nb = decompose_matrix_to_nonbranching(mat, precision);
    for (auto const &m : mats_nb) {
      if constexpr (isreal<coeff_t>()) {
        if (m.isreal()) {
          ops.push_back(NonBranchingOp<bit_t, coeff_t>(sites, m.as<arma::mat>(),
                                                       d, precision));
        } else {
          XDIAG_THROW(
              "Cannot create a real NonBranchingOp from a complex matrix.")
        }
      } else {
        ops.push_back(NonBranchingOp<bit_t, coeff_t>(
            sites, m.as<arma::cx_mat>(), d, precision));
      }
    }
  } else {
    XDIAG_THROW(
        "Cannot convert Op to NonBranchingOps. Op has no matrix defined.");
  }
  return ops;
} XDIAG_CATCH

#define INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BIT_TYPE, NUMBER_TYPE)        \
  template class NonBranchingOp<BIT_TYPE, NUMBER_TYPE>;                        \
  template std::vector<NonBranchingOp<BIT_TYPE, NUMBER_TYPE>>                  \
  non_branching_ops<BIT_TYPE, NUMBER_TYPE>(Coeff const &, Op const &,          \
                                           int64_t, double);

using namespace bits;

INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint32_t, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint64_t, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetDynamic, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic2, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic4, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic8, double);

INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint32_t, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint64_t, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetDynamic, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic2, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic4, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic8, complex);

} // namespace xdiag::matrices
