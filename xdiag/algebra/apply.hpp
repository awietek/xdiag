#pragma once

#include <xdiag/common.hpp>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#endif

namespace xdiag {

XDIAG_API void apply(OpSum const &ops, State const &v, State &w,
                     double precision = 1e-12);

template <typename coeff_t>
void apply(OpSum const &op, Spinhalf const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
           arma::Col<coeff_t> &vec_out, double precision = 1e-12);
template <typename coeff_t>
void apply(OpSum const &op, Spinhalf const &block_in,
           arma::Mat<coeff_t> const &mat_in, Spinhalf const &block_out,
           arma::Mat<coeff_t> &mat_out, double precision = 1e-12);

template <typename coeff_t>
void apply(OpSum const &op, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out, double precision = 1e-12);
template <typename coeff_t>
void apply(OpSum const &op, tJ const &block_in,
           arma::Mat<coeff_t> const &mat_in, tJ const &block_out,
           arma::Mat<coeff_t> &mat_out, double precision = 1e-12);

template <typename coeff_t>
void apply(OpSum const &op, Electron const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron const &block_out,
           arma::Col<coeff_t> &vec_out, double precision = 1e-12);
template <typename coeff_t>
void apply(OpSum const &op, Electron const &block_in,
           arma::Mat<coeff_t> const &mat_in, Electron const &block_out,
           arma::Mat<coeff_t> &mat_out, double precision = 1e-12);

#ifdef XDIAG_USE_MPI
template <typename coeff_t>
void apply(OpSum const &ops, SpinhalfDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in,
           SpinhalfDistributed const &block_out, arma::Col<coeff_t> &vec_out,
           double precision = 1e-12);

template <typename coeff_t>
void apply(OpSum const &ops, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out, double precision = 1e-12);
#endif

template <typename coeff_t>
void apply(OpSum const &ops, Block const &block_in,
           arma::Col<coeff_t> const &vec_in, Block const &block_out,
           arma::Col<coeff_t> &vec_out, double precision = 1e-12);

template <typename coeff_t>
void apply(OpSum const &ops, Block const &block_in,
           arma::Mat<coeff_t> const &mat_in, Block const &block_out,
           arma::Mat<coeff_t> &mat_out, double precision = 1e-12);

} // namespace xdiag
