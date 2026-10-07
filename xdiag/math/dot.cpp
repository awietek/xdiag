// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "dot.hpp"

#include <algorithm>
#include <vector>

#include <xdiag/utils/error.hpp>

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/mpi/allreduce.hpp>
#include <xdiag/mpi/cdot_distributed.hpp>
#endif

namespace xdiag::math {

// Chunk length of dot_chunked. A compile-time constant independent of the
// thread count and the vector length, so that chunk boundaries, the chunk sums
// and their in-order reduction are identical for any number of threads.
constexpr int64_t dot_chunk_size = 4096;

// Number of independent accumulators in the chunk loop. Breaks the dependency
// chain of a single accumulator (so the loop vectorizes without reassociation
// flags) and fixes the summation pattern in the source.
constexpr int64_t dot_nacc = 8;

template <typename coeff_t>
coeff_t dot_chunked(arma::Col<coeff_t> const &v,
                    arma::Col<coeff_t> const &w) try {
  if (v.n_rows != w.n_rows) {
    XDIAG_THROW("vector size does not match");
  }
  int64_t n = v.n_rows;
  int64_t nchunks = (n + dot_chunk_size - 1) / dot_chunk_size;
  coeff_t const *vp = v.memptr();
  coeff_t const *wp = w.memptr();
  std::vector<coeff_t> chunk_sums(nchunks, coeff_t(0));

  // Chunk sums: which thread computes a chunk does not affect its value.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int64_t c = 0; c < nchunks; ++c) {
    int64_t begin = c * dot_chunk_size;
    int64_t end = std::min(begin + dot_chunk_size, n);
    coeff_t acc[dot_nacc] = {};
    int64_t i = begin;
    for (; i + dot_nacc <= end; i += dot_nacc) {
      for (int64_t k = 0; k < dot_nacc; ++k) {
        acc[k] += conj(vp[i + k]) * wp[i + k];
      }
    }
    for (; i < end; ++i) {
      acc[0] += conj(vp[i]) * wp[i];
    }
    coeff_t s = 0;
    for (int64_t k = 0; k < dot_nacc; ++k) {
      s += acc[k];
    }
    chunk_sums[c] = s;
  }

  // In-order reduction of the chunk sums (serial, never in arrival order).
  coeff_t total = 0;
  for (int64_t c = 0; c < nchunks; ++c) {
    total += chunk_sums[c];
  }
  return total;
}
XDIAG_CATCH

template double dot_chunked(arma::vec const &, arma::vec const &);
template complex dot_chunked(arma::cx_vec const &, arma::cx_vec const &);

double dot(Block const &block, arma::vec const &v, arma::vec const &w) try {
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return dot_chunked(v, w);
#ifdef XDIAG_DISTRIBUTED
  }
#endif
}
XDIAG_CATCH

complex dot(Block const &block, arma::cx_vec const &v,
            arma::cx_vec const &w) try {
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return dot_chunked(v, w);
#ifdef XDIAG_DISTRIBUTED
  }
#endif
}
XDIAG_CATCH

template <typename coeff_t>
arma::Mat<coeff_t> matrix_dot(Block const &block, arma::Mat<coeff_t> const &V,
                              arma::Mat<coeff_t> const &W) try {
  if (V.n_rows != W.n_rows) {
    XDIAG_THROW("Input matrices do not have the same number of rows");
  }
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    int64_t L = V.n_rows;
    int64_t m = V.n_cols;
    int64_t n = W.n_cols;
    arma::Mat<coeff_t> result(m, n, arma::fill::zeros);
    for (int64_t i = 0; i < m; ++i) {
      for (int64_t j = 0; j < n; ++j) {
        arma::Col<coeff_t> cv(const_cast<coeff_t *>(V.colptr(i)), L, false,
                              true);
        arma::Col<coeff_t> cw(const_cast<coeff_t *>(W.colptr(j)), L, false,
                              true);
        result(i, j) = cdot_distributed(cv, cw);
      }
    }
    return result;
  } else {
#else
  (void)block;
#endif
    return V.t() * W;
#ifdef XDIAG_DISTRIBUTED
  }
#endif
}
XDIAG_CATCH

template arma::mat matrix_dot(Block const &, arma::mat const &,
                              arma::mat const &);
template arma::cx_mat matrix_dot(Block const &, arma::cx_mat const &,
                                 arma::cx_mat const &);

} // namespace xdiag::math
