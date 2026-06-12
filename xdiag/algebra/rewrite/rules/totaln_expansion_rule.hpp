// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// Rule to expand the site-free total particle number operator into a sum of
// local number operators over all `nsites` sites of the algebra:
//
//   TotalN -> sum_i N{i}
//
OpRule totaln_expansion_rule(int64_t nsites);

} // namespace xdiag::algebra
