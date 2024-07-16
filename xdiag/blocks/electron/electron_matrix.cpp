#include "electron_matrix.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/electron/compile.hpp>
#include <xdiag/blocks/electron/dispatch.hpp>
#include <xdiag/operators/compiler.hpp>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(OpSum const &ops, Electron const &block_in,
                              Electron const &block_out,
                              double zero_precision) try {
  int64_t m = block_out.size();
  int64_t n = block_in.dim();
  arma::Mat<coeff_t> mat(m, n, arma::fill::zeros);
  matrix_gen(mat.memptr(), ops, block_in, block_out, zero_precision);
  return mat;

} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

arma::mat matrix(OpSum const &ops, Electron const &block_in,
                 Electron const &block_out, double zero_precision) try {
  return matrix_gen<double>(ops, block_in, block_out, zero_precision);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

arma::cx_mat matrixC(OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, double zero_precision) try {
  return matrix_gen<complex>(ops, block_in, block_out, zero_precision);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template <typename coeff_t>
void matrix_gen(coeff_t *mat, OpSum const &ops, Electron const &block_in,
                Electron const &block_out, double zero_precision) try {
  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);

  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_matrix(mat, idx_in, idx_out, m, val);
  };

  electron::dispatch<coeff_t>(ops, block_in, block_out, fill, zero_precision);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

void matrix(double *mat, OpSum const &ops, Electron const &block_in,
            Electron const &block_out, double zero_precision) try {
  return matrix_gen<double>(mat, ops, block_in, block_out, zero_precision);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

void matrixC(complex *mat, OpSum const &ops, Electron const &block_in,
             Electron const &block_out, double zero_precision) try {
  return matrix_gen<complex>(mat, ops, block_in, block_out, zero_precision);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag
