#include "matrix.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/basis/electron/apply/dispatch_matrix.hpp>
#include <xdiag/basis/spinhalf/apply/dispatch_matrix.hpp>
#include <xdiag/basis/tj/apply/dispatch_matrix.hpp>
#include <xdiag/operators/compiler.hpp>

namespace xdiag {

template <typename coeff_t, class block_t>
arma::Mat<coeff_t> matrix_gen(OpSum const &ops, block_t const &block_in,
                              block_t const &block_out, double precision) try {
  int64_t m = block_out.size();
  int64_t n = block_in.dim();
  arma::Mat<coeff_t> mat(m, n, arma::fill::zeros);
  matrix(mat.memptr(), ops, block_in, block_out, precision);
  return mat;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <class block_t>
arma::mat matrix(OpSum const &ops, block_t const &block_in,
                 block_t const &block_out, double precision) try {
  return matrix_gen<double>(ops, block_in, block_out, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::mat matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                          double);
template arma::mat matrix(OpSum const &, tJ const &, tJ const &, double);
template arma::mat matrix(OpSum const &, Electron const &, Electron const &,
                          double);

template <class block_t>
arma::mat matrix(Op const &op, block_t const &block_in,
                 block_t const &block_out, double precision) try {
  OpSum ops({op});
  return matrix(ops, block_in, block_out, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::mat matrix(Op const &, Spinhalf const &, Spinhalf const &,
                          double);
template arma::mat matrix(Op const &, tJ const &, tJ const &, double);
template arma::mat matrix(Op const &, Electron const &, Electron const &,
                          double);

template <class block_t>
arma::mat matrix(OpSum const &ops, block_t const &block, double precision) try {
  return matrix(ops, block, block, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::mat matrix(OpSum const &, Spinhalf const &, double);
template arma::mat matrix(OpSum const &, tJ const &, double);
template arma::mat matrix(OpSum const &, Electron const &, double);

template <class block_t>
arma::mat matrix(Op const &op, block_t const &block, double precision) try {
  return matrix(op, block, block, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::mat matrix(Op const &, Spinhalf const &, double);
template arma::mat matrix(Op const &, tJ const &, double);
template arma::mat matrix(Op const &, Electron const &, double);

template <class block_t>
arma::cx_mat matrixC(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out, double precision) try {
  return matrix_gen<complex>(ops, block_in, block_out, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::cx_mat matrixC(OpSum const &, Spinhalf const &, Spinhalf const &,
                              double);
template arma::cx_mat matrixC(OpSum const &, tJ const &, tJ const &, double);
template arma::cx_mat matrixC(OpSum const &, Electron const &, Electron const &,
                              double);

template <class block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block_in,
                     block_t const &block_out, double precision) try {
  OpSum ops({op});
  return matrixC(ops, block_in, block_out, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::cx_mat matrixC(Op const &, Spinhalf const &, Spinhalf const &,
                              double);
template arma::cx_mat matrixC(Op const &, tJ const &, tJ const &, double);
template arma::cx_mat matrixC(Op const &, Electron const &, Electron const &,
                              double);

template <class block_t>
arma::cx_mat matrixC(OpSum const &ops, block_t const &block,
                     double precision) try {
  return matrixC(ops, block, block, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::cx_mat matrixC(OpSum const &, Spinhalf const &, double);
template arma::cx_mat matrixC(OpSum const &, tJ const &, double);
template arma::cx_mat matrixC(OpSum const &, Electron const &, double);

template <class block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block, double precision) try {
  return matrixC(op, block, block, precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::cx_mat matrixC(Op const &, Spinhalf const &, double);
template arma::cx_mat matrixC(Op const &, tJ const &, double);
template arma::cx_mat matrixC(Op const &, Electron const &, double);

// developer methods
template <typename coeff_t>
void matrix(coeff_t *mat, OpSum const &ops, Spinhalf const &block_in,
            Spinhalf const &block_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_spinhalf(ops, n_sites, precision);
  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);
  basis::spinhalf::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void matrix(double *mat, OpSum const &ops, Spinhalf const &block_in,
                     Spinhalf const &block_out, double precision);
template void matrix(complex *mat, OpSum const &ops, Spinhalf const &block_in,
                     Spinhalf const &block_out, double precision);

template <typename coeff_t>
void matrix(coeff_t *mat, OpSum const &ops, tJ const &block_in,
            tJ const &block_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_tj(ops, n_sites, precision);
  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);
  basis::tj::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void matrix(double *mat, OpSum const &ops, tJ const &block_in,
                     tJ const &block_out, double precision);
template void matrix(complex *mat, OpSum const &ops, tJ const &block_in,
                     tJ const &block_out, double precision);

template <typename coeff_t>
void matrix(coeff_t *mat, OpSum const &ops, Electron const &block_in,
            Electron const &block_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_electron(ops, n_sites, precision);
  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);
  basis::electron::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void matrix(double *mat, OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, double precision);
template void matrix(complex *mat, OpSum const &ops, Electron const &block_in,
                     Electron const &block_out, double precision);

} // namespace xdiag
