// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/armadillo.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Correlation matrix of two single-site operators:
//
//   C(i, j) = <state| Op(type1, i) Op(type2, j) |state>,   i,j = 0 .. nsites-1.
//
// `state` is assumed normalized (raw matrix elements are returned, mirroring
// ITensor's correlation_matrix). The product Op(type1,i)*Op(type2,j) must keep
// the state within its block (e.g. type1="Adag", type2="A" gives the
// single-particle density matrix); a product that changes the block's conserved
// quantum numbers is reported as a block mismatch. On a symmetry-adapted block
// the operator is symmetrized with the trivial irrep so it acts within the
// block.
//
// `correlation_matrix` returns a real matrix and throws if the result is
// complex; use `correlation_matrixC` for a complex result.
XDIAG_API arma::mat correlation_matrix(State const &state, std::string type1,
                                       std::string type2);
XDIAG_API arma::cx_mat correlation_matrixC(State const &state,
                                           std::string type1,
                                           std::string type2);

} // namespace xdiag
