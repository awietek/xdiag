#pragma once
#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/operators/bond.hpp>
#include <xdiag/operators/bondlist.hpp>

#include <xdiag/blocks/electron/electron_matrix.hpp>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.hpp>
#include <xdiag/blocks/tj/tj_matrix.hpp>

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
