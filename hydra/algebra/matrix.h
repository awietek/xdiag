#pragma once
#include "extern/armadillo/armadillo"

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>

#include <hydra/blocks/electron/electron_matrix.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>
#include <hydra/blocks/tj/tj_matrix.h>

namespace hydra {

// real matrices
template <typename block_t>
arma::mat matrix(Bond const &bond, block_t const &block_in,
                 block_t const &block_out);
template <typename block_t>
arma::mat matrix(BondList const &bonds, block_t const &block);
template <typename block_t>
arma::mat matrix(Bond const &bond, block_t const &block);

// complex matrices
template <typename block_t>
arma::cx_mat matrixC(Bond const &bond, block_t const &block_in,
                     block_t const &block_out);
template <typename block_t>
arma::cx_mat matrixC(BondList const &bonds, block_t const &block);
template <typename block_t>
arma::cx_mat matrixC(Bond const &bond, block_t const &block);

} // namespace hydra
