#include "matrix.h"

#include <hydra/blocks/electron/electron_matrix.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>
#include <hydra/blocks/tj/tj_matrix.h>

namespace hydra {

arma::Mat<double> matrix_real(BondList const &bonds, Block const &block_in,
                              Block const &block_out) {
  return std::visit(variant::overloaded{
      [&](BondList const &bonds, Spinhalf const &block_in,
          Spinhalf const &block_out) -> arma::Mat<double> {
        return matrix_real(bonds, block_in, block_out);
      },
      [&](BondList const &bonds, tJ const &block_in,
          tJ const &block_out) -> arma::Mat<double> {
        return matrix_real(bonds, block_in, block_out);
      },
      [&](BondList const &bonds, Electron const &block_in,
          Electron const &block_out) -> arma::Mat<double> {
        return matrix_real(bonds, block_in, block_out);
      },
      [&](auto const &block_in, auto const &block_out) -> arma::Mat<double> {
        Log.err("Error in matrix_real: Invalid blocks (combination)!");
        return arma::mat();
      }}, block_in, block_out);
}
arma::Mat<double> matrix_real(Bond const &bond, Block const &block_in,
                              Block const &block_out) {
  BondList bonds({bond});
  return matrix_real(bonds, block_in, block_out);
}

arma::Mat<double> matrix_real(BondList const &bonds, Block const &block) {
  return matrix_real(bonds, block, block);
}

arma::Mat<double> matrix_real(Bond const &bond, Block const &block) {
  BondList bonds({bond});
  return matrix_real(bonds, block);
}

arma::Mat<complex> matrix_cplx(BondList const &bonds, Block const &block_in,
                               Block const &block_out) {
  return std::visit(variant::overloaded{
      [&](BondList const &bonds, Spinhalf const &block_in,
          Spinhalf const &block_out) -> arma::Mat<complex> {
        return matrix_cplx(bonds, block_in, block_out);
      },
      [&](BondList const &bonds, tJ const &block_in,
          tJ const &block_out) -> arma::Mat<complex> {
        return matrix_cplx(bonds, block_in, block_out);
      },
      [&](BondList const &bonds, Electron const &block_in,
          Electron const &block_out) -> arma::Mat<complex> {
        return matrix_cplx(bonds, block_in, block_out);
      },
      [&](auto const &block_in, auto const &block_out) {
        Log.err("Error in matrix_real: Invalid blocks (combination)!");
        return arma::cx_mat();
      }}, block_in, block_out);
}
arma::Mat<complex> matrix_cplx(Bond const &bond, Block const &block_in,
                               Block const &block_out) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block_in, block_out);
}

arma::Mat<complex> matrix_cplx(BondList const &bonds, Block const &block) {
  return matrix_cplx(bonds, block, block);
}

arma::Mat<complex> matrix_cplx(Bond const &bond, Block const &block) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block);
}

arma::Mat<complex> matrix(BondList const &bonds, Block const &block_in,
                          Block const &block_out) {
  return matrix_cplx(bonds, block_in, block_out);
}

arma::Mat<complex> matrix(Bond const &bond, Block const &block_in,
                          Block const &block_out) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block_in, block_out);
}

arma::Mat<complex> matrix(BondList const &bonds, Block const &block) {
  return matrix_cplx(bonds, block, block);
}

arma::Mat<complex> matrix(Bond const &bond, Block const &block) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block);
}

} // namespace hydra
