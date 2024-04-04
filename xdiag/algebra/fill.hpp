#pragma once

#include <xdiag/common.hpp>

namespace xdiag {

template <typename coeff_t>
constexpr inline void fill_matrix(coeff_t *memptr, int64_t idx_in,
                                  int64_t idx_out, int64_t m, coeff_t val) {
  memptr[idx_out + idx_in * m] += val;
}

template <typename coeff_t>
constexpr inline void fill_matrix(arma::Mat<coeff_t> &mat, int64_t idx_in,
                                  int64_t idx_out, coeff_t val) {
  fill_matrix(mat.memptr(), idx_in, idx_out, mat.n_rows, val);
}

template <typename coeff_t>
inline void fill_apply(coeff_t const *vec_in, coeff_t *vec_out, int64_t idx_in,
                       int64_t idx_out, coeff_t val) {
  // Atomic update to avoid multiple threads writing to the same address
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
inline void fill_apply(arma::Col<coeff_t> const &vec_in,
                       arma::Col<coeff_t> &vec_out, int64_t idx_in,
                       int64_t idx_out, coeff_t val) {
  fill_apply(vec_in.memptr(), vec_out.memptr(), idx_in, idx_out, val);
}

} // namespace xdiag
