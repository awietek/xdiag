#include "tj_matrix.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/tj/compile.hpp>
#include <xdiag/blocks/tj/dispatch.hpp>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, tJ const &block_in,
                              tJ const &block_out, double zero_precision) try {
  int64_t m = block_out.size();
  int64_t n = block_in.dim();
  arma::Mat<coeff_t> mat(m, n, arma::fill::zeros);
  matrix_gen(mat.memptr(), bonds, block_in, block_out, zero_precision);
  return mat;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

arma::mat matrix(BondList const &bonds, tJ const &block_in, tJ const &block_out,
                 double zero_precision) try {
  return matrix_gen<double>(bonds, block_in, block_out, zero_precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

arma::cx_mat matrixC(BondList const &bonds, tJ const &block_in,
                     tJ const &block_out, double zero_precision) try {
  return matrix_gen<complex>(bonds, block_in, block_out, zero_precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t>
void matrix_gen(coeff_t *mat, BondList const &bonds, tJ const &block_in,
                tJ const &block_out, double zero_precision) try {

  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);

  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_matrix(mat, idx_in, idx_out, m, val);
  };

  tj::dispatch<coeff_t>(bonds, block_in, block_out, fill, zero_precision);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void matrix(double *mat, BondList const &bonds, tJ const &block_in,
            tJ const &block_out, double zero_precision) try {
  return matrix_gen<double>(mat, bonds, block_in, block_out, zero_precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void matrixC(complex *mat, BondList const &bonds, tJ const &block_in,
             tJ const &block_out, double zero_precision) try {
  return matrix_gen<complex>(mat, bonds, block_in, block_out, zero_precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
