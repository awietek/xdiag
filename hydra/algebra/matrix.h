#pragma once
#include "extern/armadillo/armadillo"

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

arma::Mat<double> matrix_real(BondList const &bond, Block const &block_in,
                              Block const &block_out);
arma::Mat<double> matrix_real(Bond const &bond, Block const &block_in,
                              Block const &block_out);
arma::Mat<double> matrix_real(BondList const &bonds, Block const &block);
arma::Mat<double> matrix_real(Bond const &bond, Block const &block);

arma::Mat<complex> matrix_cplx(BondList const &bond, Block const &block_in,
                               Block const &block_out);
arma::Mat<complex> matrix_cplx(Bond const &bond, Block const &block_in,
                               Block const &block_out);
arma::Mat<complex> matrix_cplx(BondList const &bonds, Block const &block);
arma::Mat<complex> matrix_cplx(Bond const &bond, Block const &block);

arma::Mat<complex> matrix(BondList const &bonds, Block const &block_in,
                          Block const &block_out);
arma::Mat<complex> matrix(Bond const &bond, Block const &block_in,
                          Block const &block_out);
arma::Mat<complex> matrix(BondList const &bonds, Block const &block);
arma::Mat<complex> matrix(Bond const &bond, Block const &block);
  
} // namespace hydra
