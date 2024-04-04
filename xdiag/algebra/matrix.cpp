#include "matrix.hpp"

#include <xdiag/blocks/electron/electron_matrix.hpp>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.hpp>
#include <xdiag/blocks/tj/tj_matrix.hpp>

namespace xdiag {

template <typename block_t>
arma::mat matrix(Bond const &bond, block_t const &block_in,
                 block_t const &block_out) try {
  BondList bonds({bond});
  return matrix(bonds, block_in, block_out);
} catch (...) {
  XDiagRethrow("Unable to create real matrix");
  return arma::mat();
}

template arma::mat matrix(Bond const &, Spinhalf const &, Spinhalf const &);
template arma::mat matrix(Bond const &, tJ const &, tJ const &);
template arma::mat matrix(Bond const &, Electron const &, Electron const &);

template <typename block_t>
arma::mat matrix(BondList const &bonds, block_t const &block) try {
  return matrix(bonds, block, block);
} catch (...) {
  XDiagRethrow("Unable to create real matrix");
  return arma::mat();
}
template arma::mat matrix(BondList const &, Spinhalf const &);
template arma::mat matrix(BondList const &, tJ const &);
template arma::mat matrix(BondList const &, Electron const &);

template <typename block_t>
arma::mat matrix(Bond const &bond, block_t const &block) try {
  return matrix(bond, block, block);
} catch (...) {
  XDiagRethrow("Unable to create real matrix");
  return arma::mat();
}
template arma::mat matrix(Bond const &, Spinhalf const &);
template arma::mat matrix(Bond const &, tJ const &);
template arma::mat matrix(Bond const &, Electron const &);

template <typename block_t>
arma::cx_mat matrixC(Bond const &bond, block_t const &block_in,
                     block_t const &block_out) try {
  BondList bonds({bond});
  return matrixC(bonds, block_in, block_out);
} catch (...) {
  XDiagRethrow("Unable to create complex matrix");
  return arma::cx_mat();
}
template arma::cx_mat matrixC(Bond const &, Spinhalf const &, Spinhalf const &);
template arma::cx_mat matrixC(Bond const &, tJ const &, tJ const &);
template arma::cx_mat matrixC(Bond const &, Electron const &, Electron const &);

template <typename block_t>
arma::cx_mat matrixC(BondList const &bonds, block_t const &block) try {
  return matrixC(bonds, block, block);
} catch (...) {
  XDiagRethrow("Unable to create complex matrix");
  return arma::cx_mat();
}
template arma::cx_mat matrixC(BondList const &, Spinhalf const &);
template arma::cx_mat matrixC(BondList const &, tJ const &);
template arma::cx_mat matrixC(BondList const &, Electron const &);

template <typename block_t>
arma::cx_mat matrixC(Bond const &bond, block_t const &block) try {
  return matrixC(bond, block, block);
} catch (...) {
  XDiagRethrow("Unable to create complex matrix");
  return arma::cx_mat();
}
template arma::cx_mat matrixC(Bond const &, Spinhalf const &);
template arma::cx_mat matrixC(Bond const &, tJ const &);
template arma::cx_mat matrixC(Bond const &, Electron const &);

} // namespace xdiag
