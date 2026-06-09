// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// Bosonic expansion rules used by the matrix algebra. They expand the compound
// particle operators into products of the elementary ladder/number operators
// (Adag, A, N), which are subsequently converted to explicit matrices:
//
//   Hop{i,j} -> -Adag{i}*A{j} - Adag{j}*A{i}  (minus sign, field convention,
//                                              matching the fermionic Hop)
//   HubbardU -> (1/2) sum_i N{i}*(N{i} - 1)   (Bose-Hubbard on-site interaction)
//
// HubbardU has no sites and sums over all `nsites` sites of the algebra, so
// `nsites` must be supplied here.
std::vector<OpRule> boson_expansion_rules(int64_t nsites);

} // namespace xdiag::algebra
