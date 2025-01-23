#pragma once
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

namespace xdiag {

template <class block_t>
XDIAG_API arma::mat matrix(OpSum const &ops, block_t const &block);
template <class block_t>
XDIAG_API arma::mat matrix(Op const &op, block_t const &block);

template <class block_t>
XDIAG_API arma::mat matrix(OpSum const &ops, block_t const &block_in,
                           block_t const &block_out);
template <class block_t>
XDIAG_API arma::mat matrix(Op const &op, block_t const &block_in,
                           block_t const &block_out);

template <class block_t>
XDIAG_API arma::cx_mat matrixC(OpSum const &ops, block_t const &block);
template <class block_t>
XDIAG_API arma::cx_mat matrixC(Op const &op, block_t const &block);

template <class block_t>
XDIAG_API arma::cx_mat matrixC(OpSum const &ops, block_t const &block_in,
                               block_t const &block_out);
template <class block_t>
XDIAG_API arma::cx_mat matrixC(Op const &op, block_t const &block_in,
                               block_t const &block_out);

// developer methods
template <typename coeff_t, class block_t>
XDIAG_API void matrix(coeff_t *mat, OpSum const &ops, block_t const &block_in,
                      block_t const &block_out);
template <typename coeff_t, class block_t>
XDIAG_API void matrix(coeff_t *mat, Op const &ops, block_t const &block_in,
                      block_t const &block_out);

template <typename coeff_t, class block_t>
XDIAG_API void matrix(coeff_t *mat, OpSum const &ops, block_t const &block);
template <typename coeff_t, class block_t>
XDIAG_API void matrix(coeff_t *mat, Op const &ops, block_t const &block);

} // namespace xdiag
