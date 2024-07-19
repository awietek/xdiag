#include "matrix.hpp"

namespace xdiag {

template <typename block_t>
arma::mat matrix(Op const &op, block_t const &block_in,
                 block_t const &block_out) try {
  OpSum ops({op});
  return matrix(ops, block_in, block_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return arma::mat();
}

template arma::mat matrix(Op const &, Spinhalf const &, Spinhalf const &);
template arma::mat matrix(Op const &, tJ const &, tJ const &);
template arma::mat matrix(Op const &, Electron const &, Electron const &);

template <typename block_t>
arma::mat matrix(OpSum const &ops, block_t const &block) try {
  return matrix(ops, block, block);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return arma::mat();
}
template arma::mat matrix(OpSum const &, Spinhalf const &);
template arma::mat matrix(OpSum const &, tJ const &);
template arma::mat matrix(OpSum const &, Electron const &);

template <typename block_t>
arma::mat matrix(Op const &op, block_t const &block) try {
  return matrix(op, block, block);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return arma::mat();
}
template arma::mat matrix(Op const &, Spinhalf const &);
template arma::mat matrix(Op const &, tJ const &);
template arma::mat matrix(Op const &, Electron const &);

template <typename block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block_in,
                     block_t const &block_out) try {
  OpSum ops({op});
  return matrixC(ops, block_in, block_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return arma::cx_mat();
}
template arma::cx_mat matrixC(Op const &, Spinhalf const &, Spinhalf const &);
template arma::cx_mat matrixC(Op const &, tJ const &, tJ const &);
template arma::cx_mat matrixC(Op const &, Electron const &, Electron const &);

template <typename block_t>
arma::cx_mat matrixC(OpSum const &ops, block_t const &block) try {
  return matrixC(ops, block, block);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return arma::cx_mat();
}
template arma::cx_mat matrixC(OpSum const &, Spinhalf const &);
template arma::cx_mat matrixC(OpSum const &, tJ const &);
template arma::cx_mat matrixC(OpSum const &, Electron const &);

template <typename block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block) try {
  return matrixC(op, block, block);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return arma::cx_mat();
}
template arma::cx_mat matrixC(Op const &, Spinhalf const &);
template arma::cx_mat matrixC(Op const &, tJ const &);
template arma::cx_mat matrixC(Op const &, Electron const &);

} // namespace xdiag
