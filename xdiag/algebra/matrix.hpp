#pragma once
#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

#include <xdiag/blocks/electron/matrix.hpp>
#include <xdiag/blocks/spinhalf/matrix.hpp>
#include <xdiag/blocks/tj/matrix.hpp>

namespace xdiag {

// real matrices
template <typename block_t>
arma::mat matrix(Op const &op, block_t const &block_in,
                 block_t const &block_out);
template <typename block_t>
arma::mat matrix(OpSum const &ops, block_t const &block);
template <typename block_t>
arma::mat matrix(Op const &op, block_t const &block);

// complex matrices
template <typename block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block_in,
                     block_t const &block_out);
template <typename block_t>
arma::cx_mat matrixC(OpSum const &ops, block_t const &block);
template <typename block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block);

} // namespace xdiag
