// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/distributed/electron_distributed.hpp>
#include <xdiag/blocks/distributed/spinhalf_distributed.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/opsum.hpp>

// Matrix-free, MPI-aware apply for the distributed blocks. One overload per
// distributed block type; each is defined and explicitly instantiated in that
// block's kernels.cpp (which pulls in the heavy MPI term kernels). apply.cpp
// only needs these declarations.
namespace xdiag::kernels {

template <typename coeff_t>
void apply_distributed(OpSum const &ops, SpinhalfDistributed const &block_in,
                       arma::Col<coeff_t> const &vec_in,
                       SpinhalfDistributed const &block_out,
                       arma::Col<coeff_t> &vec_out);

template <typename coeff_t>
void apply_distributed(OpSum const &ops, tJDistributed const &block_in,
                       arma::Col<coeff_t> const &vec_in,
                       tJDistributed const &block_out,
                       arma::Col<coeff_t> &vec_out);

template <typename coeff_t>
void apply_distributed(OpSum const &ops, ElectronDistributed const &block_in,
                       arma::Col<coeff_t> const &vec_in,
                       ElectronDistributed const &block_out,
                       arma::Col<coeff_t> &vec_out);

template <typename coeff_t>
void apply_distributed(OpSum const &ops, Block const &block_in,
                       arma::Mat<coeff_t> const &vec_in, Block const &block_out,
                       arma::Mat<coeff_t> &vec_out);

} // namespace xdiag::kernels
#endif
