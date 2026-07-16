// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <utility>
#include <vector>

#include <tests/catch.hpp>

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/correlation_matrix.hpp>
#include <xdiag/states/expect.hpp>

namespace xdiag {

// Compute the ground state on a distributed block and its non-distributed
// counterpart for the same (Hermitian, unique-ground-state) OpSum, and check
// that expect / correlation_matrix observables agree. These are phase-invariant
// (<psi|O|psi>), so they match whenever the ground state is non-degenerate --
// shared by the tJ / electron / spinhalf distributed test files.
inline void compare_observables(
    OpSum const &ops, Block block, Block block_dist,
    std::vector<std::string> const &onesite,
    std::vector<std::pair<std::string, std::string>> const &twosite) {
  auto [e0, psi] = eig0(ops, block);
  auto [e0d, psid] = eig0(ops, block_dist);
  REQUIRE(std::abs(e0 - e0d) < 1e-8);
  for (auto const &type : onesite) {
    arma::cx_vec x = expectC(psi, type);
    arma::cx_vec xd = expectC(psid, type);
    REQUIRE(arma::norm(x - xd) < 1e-6);
  }
  for (auto const &[t1, t2] : twosite) {
    arma::cx_mat c = correlation_matrixC(psi, t1, t2);
    arma::cx_mat cd = correlation_matrixC(psid, t1, t2);
    REQUIRE(arma::norm(c - cd) < 1e-6);
  }
}

} // namespace xdiag
