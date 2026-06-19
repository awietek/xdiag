// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/lambda_details.hpp>

namespace xdiag::matrices {

// ---------------------------------------------------------------------------
// Dense matrix fill (column-major, no atomics needed: each (idx_in,idx_out)
// maps to a unique element so concurrent writes cannot collide).
// ---------------------------------------------------------------------------

template <typename coeff_t>
constexpr void fill_matrix(coeff_t *mat, int64_t m, int64_t idx_in,
                           int64_t idx_out, coeff_t val) {
  mat[idx_out + idx_in * m] += val; // column-major order
};

// ---------------------------------------------------------------------------
// Apply fill  (vec_out[idx_out] += val * vec_in[idx_in]).
// Atomics required because multiple threads may contribute to the same row.
// ---------------------------------------------------------------------------

template <typename coeff_t>
inline void fill_apply(coeff_t const *vec_in, coeff_t *vec_out, int64_t idx_in,
                       int64_t idx_out, coeff_t val) {
#ifdef _OPENMP
  if constexpr (isreal<coeff_t>()) {
    coeff_t x = val * vec_in[idx_in];
#pragma omp atomic update
    vec_out[idx_out] += x;
  } else {
    complex x = val * vec_in[idx_in];
    double *r = &reinterpret_cast<double(&)[2]>(vec_out[idx_out])[0];
    double *i = &reinterpret_cast<double(&)[2]>(vec_out[idx_out])[1];
#pragma omp atomic update
    *r += x.real();
#pragma omp atomic update
    *i += x.imag();
  }
#else
  vec_out[idx_out] += val * vec_in[idx_in];
#endif
}

template <typename coeff_t>
constexpr void fill_apply(arma::Col<coeff_t> const &vec_in,
                          arma::Col<coeff_t> &vec_out, int64_t idx_in,
                          int64_t idx_out, coeff_t val) {
  fill_apply(vec_in.memptr(), vec_out.memptr(), idx_in, idx_out, val);
}

template <typename coeff_t>
constexpr void fill_apply(arma::Mat<coeff_t> const &mat_in,
                          arma::Mat<coeff_t> &mat_out, int64_t idx_in,
                          int64_t idx_out, coeff_t val) {
  for (int i = 0; i < mat_in.n_cols; i++) {
    fill_apply(mat_in.colptr(i), mat_out.colptr(i), idx_in, idx_out, val);
  }
}

// ---------------------------------------------------------------------------
// COO NNZ counting.
// Two overloads: the compiler (via #ifdef _OPENMP in the caller) selects the
// right one — per-thread counting (OMP) or a single total counter (serial).
// ---------------------------------------------------------------------------

// OMP: count into per-thread slot identified by num_thread.
#ifdef _OPENMP
inline void fill_coo_count(std::vector<int64_t> &nnz_thread, int num_thread) {
  ++nnz_thread[num_thread];
}
#else
// Serial: count into a single accumulator.
inline void fill_coo_count(int64_t &nnz) { ++nnz; }
#endif

// ---------------------------------------------------------------------------
// COO fill.
// Two overloads: OMP uses a per-thread offset array; serial uses a single
// sequential counter k.
// ---------------------------------------------------------------------------

// OMP: each thread writes to its own pre-computed offset region.
#ifdef _OPENMP
template <typename idx_t, typename coeff_t>
inline void fill_coo(idx_t *rows, idx_t *cols, coeff_t *data,
                     std::vector<int64_t> &offset, int64_t idx_in,
                     int64_t idx_out, coeff_t val, idx_t i0, int num_thread) {
  int64_t k = offset[num_thread]++;
  rows[k] = (idx_t)(idx_out + i0);
  cols[k] = (idx_t)(idx_in + i0);
  data[k] = val;
}
#else
// Serial: sequential slot assignment via incrementing counter k.
template <typename idx_t, typename coeff_t>
inline void fill_coo(idx_t *rows, idx_t *cols, coeff_t *data, int64_t &k,
                     int64_t idx_in, int64_t idx_out, coeff_t val, idx_t i0) {
  rows[k] = (idx_t)(idx_out + i0);
  cols[k] = (idx_t)(idx_in + i0);
  data[k] = val;
  ++k;
}
#endif

// ---------------------------------------------------------------------------
// CSR NNZ counting (per-row).
// Uniform across OMP/serial: omp atomic update handles concurrent increments.
// ---------------------------------------------------------------------------

inline void fill_csr_count(std::vector<int64_t> &n_elements_in_row,
                           int64_t idx_out) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
  ++n_elements_in_row[idx_out];
}

// ---------------------------------------------------------------------------
// CSR fill.
// Uniform across OMP/serial: omp atomic capture claims a unique slot within
// the row; the offset is shared across all threads writing to the same row.
// ---------------------------------------------------------------------------

template <typename idx_t, typename coeff_t>
inline void fill_csr(std::vector<int64_t> &offset, idx_t *col, coeff_t *data,
                     int64_t idx_in, int64_t idx_out, coeff_t val, idx_t i0) {
  int64_t k;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
  k = offset[idx_out]++;
  col[k] = (idx_t)(idx_in + i0);
  data[k] = val;
}

} // namespace xdiag::matrices

// ---------------------------------------------------------------------------
// XDIAG_FILL macro: dispatches fill(idx_in, idx_out, coeff) or
// fill(idx_in, idx_out, coeff, num_thread) depending on lambda arity.
// Requires "fill" and (under OpenMP) "num_thread" to be in scope.
// ---------------------------------------------------------------------------
#ifdef _OPENMP
#define XDIAG_FILL(idx_in, idx_out, coeff)                                     \
  if constexpr (::xdiag::utils::lambda_details<                                \
                    decltype(fill)>::argument_count == 3) {                    \
    fill(idx_in, idx_out, coeff);                                              \
  } else {                                                                     \
    fill(idx_in, idx_out, coeff, num_thread);                                  \
  }
#else
#define XDIAG_FILL(idx_in, idx_out, coeff) fill(idx_in, idx_out, coeff)
#endif
