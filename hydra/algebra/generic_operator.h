#pragma once
#include <hydra/blocks/electron/terms/apply_terms_dispatch.h>
#include <hydra/blocks/spinhalf/terms/apply_terms_dispatch.h>
#include <hydra/blocks/tj/terms/apply_terms_dispatch.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <typename coeff_t, typename block_t, class compile_f, class fill_f>
void generic_operator(BondList const &bonds, block_t const &block_in,
                      block_t const &block_out, compile_f compile,
                      fill_f fill) try {
  BondList bonds_compiled = compile(bonds, 1e-12);
  if (bonds_compiled.n_sites() > block_in.n_sites()) {
    HydraThrow(std::logic_error,
               "number of sites in BondList larger than in input block");
  }

  if ((isreal<coeff_t>()) && (bonds_compiled.iscomplex())) {
    HydraThrow(std::logic_error,
               "Trying to create a real operator from an "
               "intrisically complex BondList, consider using complex function "
               "instead (e.g. matrixC, applyC, etc)");
  }
  apply_terms_dispatch<coeff_t>(bonds_compiled, block_in, block_out, fill);
} catch (...) {
  HydraRethrow("Cannot evaluate generic operation");
}

template <typename coeff_t, typename block_t, class compile_f>
void generic_matrix(coeff_t *memptr, BondList const &bonds,
                    block_t const &block_in, block_t const &block_out,
                    compile_f compile) try {
  int64_t dim_in = block_in.size();
  int64_t dim_out = block_out.size();
  std::fill(memptr, memptr + dim_in * dim_out, 0.0);
  auto fill = [memptr, &dim_out](int64_t idx_out, int64_t idx_in, coeff_t val) {
    memptr[idx_out + idx_in * dim_out] += val;
  };
  generic_operator<coeff_t>(bonds, block_in, block_out, compile, fill);
} catch (...) {
  HydraRethrow("Cannot evaluate generic matrix");
}

template <typename coeff_t, typename block_t, class compile_f>
arma::Mat<coeff_t>
generic_matrix(BondList const &bonds, block_t const &block_in,
               block_t const &block_out, compile_f compile) try {
  int64_t dim_in = block_in.size();
  int64_t dim_out = block_out.size();
  arma::Mat<coeff_t> mat(dim_out, dim_in);
  generic_matrix<coeff_t>(mat.memptr(), bonds, block_in, block_out, compile);
  return mat;
} catch (...) {
  HydraRethrow("Cannot evaluate generic matrix");
  return arma::Mat<coeff_t>();
}

template <typename coeff_t, typename block_t, class compile_f>
void generic_apply(BondList const &bonds, block_t const &block_in,
                   coeff_t const *vec_in, block_t const &block_out,
                   coeff_t *vec_out, compile_f compile) try {

  int64_t dim_out = block_out.size();
  std::fill(vec_out, vec_out + dim_out, 0.0);

  auto fill = [&vec_out, &vec_in](int64_t idx_out, int64_t idx_in,
                                  coeff_t val) {
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
  };
  generic_operator<coeff_t>(bonds, block_in, block_out, compile, fill);
} catch (...) {
  HydraRethrow("Cannot evaluate generic apply");
}

template <typename coeff_t, typename block_t, class compile_f>
void generic_apply(BondList const &bonds, block_t const &block_in,
                   arma::Col<coeff_t> const &vec_in, block_t const &block_out,
                   arma::Col<coeff_t> &vec_out, compile_f compile) try {
  generic_apply<coeff_t>(bonds, block_in, vec_in.memptr(), block_out,
                         vec_out.memptr(), compile);
} catch (...) {
  HydraRethrow("Cannot evaluate generic apply");
}
} // namespace hydra
