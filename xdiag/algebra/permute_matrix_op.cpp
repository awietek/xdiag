// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "permute_matrix_op.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

#include <xdiag/math/ipow.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// Permutes rows and columns of a square matrix according to a site permutation,
// for a system with local dimension d per site. Matrix size must be d^nsites
// where nsites = perm.size().
template <typename T>
static arma::Mat<T> permute_matrix(arma::Mat<T> const &mat,
                                   Permutation const &perm, int64_t d) try {
  int64_t m = mat.n_rows;
  int64_t n = mat.n_cols;
  int64_t nsites = perm.size();
  int64_t expected = math::ipow(d, nsites);
  if (m != n || m != expected) {
    XDIAG_THROW(fmt::format(
        "Matrix dimensions ({}, {}) inconsistent with d={} and nsites={}: "
        "expected {}x{}.",
        m, n, d, nsites, expected, expected));
  }

  // Precompute the strides used by the inverse permutation: digit k of the
  // input index ends up multiplied by stride[k] in the output index.
  Permutation perm_inv = perm.inv();
  std::vector<int64_t> stride(nsites);
  for (int64_t i = 0; i < nsites; ++i) {
    stride[i] = math::ipow(d, perm_inv[i]);
  }

  auto permuted = [&](int64_t idx) {
    int64_t out = 0;
    for (int64_t k = 0; idx; ++k, idx /= d) {
      out += (idx % d) * stride[k];
    }
    return out;
  };

  arma::Mat<T> mat_permuted(m, n, arma::fill::zeros);
  for (int64_t i = 0; i < m; ++i) {
    int64_t ip = permuted(i);
    for (int64_t j = 0; j < m; ++j) {
      int64_t jp = permuted(j);
      mat_permuted(ip, jp) = mat(i, j);
    }
  }
  return mat_permuted;
} XDIAG_CATCH

Op permute_matrix_op(Op const &op, int64_t d) try {
  std::vector<int64_t> const &sites = op.sites();
  int64_t n = (int64_t)sites.size();

  // Compute the permutation that sorts sites in ascending order
  std::vector<int64_t> perm_vector(n);
  std::iota(perm_vector.begin(), perm_vector.end(), 0);
  std::sort(perm_vector.begin(), perm_vector.end(),
            [&](int64_t a, int64_t b) { return sites[a] < sites[b]; });
  Permutation perm(perm_vector);

  // Apply permutation to site list
  std::vector<int64_t> sites_sorted(n);
  for (int64_t i = 0; i < n; ++i) {
    sites_sorted[i] = sites[perm_vector[i]];
  }

  Matrix const &mat = op.matrix();
  if (mat.isreal()) {
    return Op("Matrix", sites_sorted,
              permute_matrix(mat.as<arma::mat>(), perm, d));
  } else {
    return Op("Matrix", sites_sorted,
              permute_matrix(mat.as<arma::cx_mat>(), perm, d));
  }
} XDIAG_CATCH

} // namespace xdiag::operators
