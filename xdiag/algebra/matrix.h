#pragma once
#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/operators/bond.h>
#include <xdiag/operators/bondlist.h>

#include <xdiag/blocks/electron/electron_matrix.h>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.h>
#include <xdiag/blocks/tj/tj_matrix.h>

namespace xdiag {

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

} // namespace xdiag
