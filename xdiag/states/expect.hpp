// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/armadillo.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Expectation value of a single-site operator `type` on every site:
//
//   result(i) = <state| Op(type, i) |state>,   i = 0 .. nsites-1.
//
// `state` is assumed normalized (the raw matrix element is returned, mirroring
// ITensor's expect). `type` must be a valid single-site Op type for the state's
// block (e.g. "N", "Sz", "Ntot") and must keep the state within its block; an
// operator that changes the block's conserved quantum numbers makes
// <state|Op|state> identically zero by symmetry, and is reported as a block
// mismatch rather than silently returning zero.
//
// `expect` returns a real vector and throws if the result is complex; use
// `expectC` for a complex result.
XDIAG_API arma::vec expect(State const &state, std::string type);
XDIAG_API arma::cx_vec expectC(State const &state, std::string type);

} // namespace xdiag
