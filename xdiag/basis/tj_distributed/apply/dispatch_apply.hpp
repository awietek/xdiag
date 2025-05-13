// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/tj_distributed.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, tJDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out);

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, tJDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out);

} // namespace xdiag::basis
