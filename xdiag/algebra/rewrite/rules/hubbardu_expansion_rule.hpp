// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// Rule to expand the site-free Hubbard interaction into a sum of local
// double-occupancy operators over all `nsites` sites of the electron block:
//
//   HubbardU -> sum_i Nupdn{i}
//
// Used by the electron symmetry algebra. The electron implementation algebra
// keeps HubbardU as a protected type backed by a dedicated popcount kernel.
OpRule hubbardu_expansion_rule(int64_t nsites);

} // namespace xdiag::algebra
