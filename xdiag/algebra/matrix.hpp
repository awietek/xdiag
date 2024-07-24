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
arma::mat matrix(OpSum const &ops, block_t const &block_in,
                 block_t const &block_out, double precision = 1e-12);
template <class block_t>
arma::mat matrix(Op const &op, block_t const &block_in,
                 block_t const &block_out, double precision = 1e-12);
template <class block_t>
arma::mat matrix(OpSum const &ops, block_t const &block,
                 double precision = 1e-12);
template <class block_t>
arma::mat matrix(Op const &op, block_t const &block, double precision = 1e-12);

template <class block_t>
arma::cx_mat matrixC(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out, double precision = 1e-12);
template <class block_t>
arma::cx_mat matrixC(Op const &ops, block_t const &block_in,
                     block_t const &block_out, double precision = 1e-12);
template <class block_t>
arma::cx_mat matrixC(OpSum const &ops, block_t const &block,
                     double precision = 1e-12);
template <class block_t>
arma::cx_mat matrixC(Op const &op, block_t const &block,
                     double precision = 1e-12);

// developer methods
template <typename coeff_t>
void matrix(coeff_t *mat, OpSum const &ops, Spinhalf const &block_in,
            Spinhalf const &block_out, double precision = 1e-12);

template <typename coeff_t>
void matrix(coeff_t *mat, OpSum const &ops, tJ const &block_in,
            tJ const &block_out, double precision = 1e-12);

template <typename coeff_t>
void matrix(coeff_t *mat, OpSum const &ops, Electron const &block_in,
            Electron const &block_out, double precision = 1e-12);

} // namespace xdiag
