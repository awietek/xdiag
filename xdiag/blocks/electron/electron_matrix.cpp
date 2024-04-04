#include "electron_matrix.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/electron/compile.hpp>
#include <xdiag/blocks/electron/dispatch.hpp>
#include <xdiag/operators/compiler.hpp>

namespace xdiag {

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, Electron const &block_in,
                              Electron const &block_out) try {
  BondList bondsc = electron::compile(bonds, 1e-12);

  int64_t m = block_out.size();
  int64_t n = block_in.dim();
  arma::Mat<coeff_t> mat(m, n, arma::fill::zeros);
  matrix_gen(mat.memptr(), bonds, block_in, block_out);
  return mat;

} catch (...) {
  XDiagRethrow("Cannot create matrix from \"Electron\" block");
  return arma::Mat<coeff_t>();
}

arma::mat matrix(BondList const &bonds, Electron const &block_in,
                 Electron const &block_out) try {
  return matrix_gen<double>(bonds, block_in, block_out);
} catch (...) {
  XDiagRethrow("Cannot create matrix from \"Electron\" block");
  return arma::mat();
}

arma::cx_mat matrixC(BondList const &bonds, Electron const &block_in,
                     Electron const &block_out) try {
  return matrix_gen<complex>(bonds, block_in, block_out);
} catch (...) {
  XDiagRethrow("Cannot create matrix from \"Electron\" block");
  return arma::cx_mat();
}

template <typename coeff_t>
void matrix_gen(coeff_t *mat, BondList const &bonds, Electron const &block_in,
                Electron const &block_out) try {
  BondList bondsc = electron::compile(bonds, 1e-12);

  int64_t m = block_out.size();
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_matrix(mat, idx_in, idx_out, m, val);
  };

  electron::dispatch<coeff_t>(bondsc, block_in, block_out, fill);
} catch (...) {
  XDiagRethrow("Cannot create matrix from \"Electron\" block");
}

void matrix(double *mat, BondList const &bonds, Electron const &block_in,
            Electron const &block_out) try {
  return matrix_gen<double>(mat, bonds, block_in, block_out);
} catch (...) {
  XDiagRethrow("Cannot create matrix from \"Electron\" block");
}

void matrixC(complex *mat, BondList const &bonds, Electron const &block_in,
             Electron const &block_out) try {
  return matrix_gen<complex>(mat, bonds, block_in, block_out);
} catch (...) {
  XDiagRethrow("Cannot create matrix from \"Electron\" block");
}

} // namespace xdiag
