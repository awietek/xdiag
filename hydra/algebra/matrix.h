#pragma once
#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <class Block>
inline arma::Mat<double> matrix_real(Bond const &bond, Block const &block_in,
                                     Block const &block_out) {
  BondList bonds({bond});
  return matrix_real(bonds, block_in, block_out);
}

template <class Block>
inline arma::Mat<double> matrix_real(BondList const &bonds,
                                     Block const &block) {
  return matrix_real(bonds, block, block);
}

template <class Block>
inline arma::Mat<double> matrix_real(Bond const &bond, Block const &block) {
  BondList bonds({bond});
  return matrix_real(bonds, block);
}

template <class Block>
inline arma::Mat<complex> matrix_cplx(Bond const &bond, Block const &block_in,
                                      Block const &block_out) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block_in, block_out);
}

template <class Block>
inline arma::Mat<complex> matrix_cplx(BondList const &bonds,
                                      Block const &block) {
  return matrix_cplx(bonds, block, block);
}

template <class Block>
inline arma::Mat<complex> matrix_cplx(Bond const &bond, Block const &block) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block);
}

template <class Block>
inline arma::Mat<complex> matrix(BondList const &bonds, Block const &block_in,
                                 Block const &block_out) {
  return matrix_cplx(bonds, block_in, block_out);
}

template <class Block>
inline arma::Mat<complex> matrix(Bond const &bond, Block const &block_in,
                                 Block const &block_out) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block_in, block_out);
}

template <class Block>
inline arma::Mat<complex> matrix(BondList const &bonds, Block const &block) {
  return matrix_cplx(bonds, block, block);
}

template <class Block>
inline arma::Mat<complex> matrix(Bond const &bond, Block const &block) {
  BondList bonds({bond});
  return matrix_cplx(bonds, block);
}

} // namespace hydra
