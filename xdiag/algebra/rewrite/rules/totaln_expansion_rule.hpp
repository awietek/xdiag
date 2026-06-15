// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// Rule to expand the site-free total particle number operator into a sum of
// local number operators over all `nsites` sites of the algebra:
//
//   TotalN -> sum_i <local_type>{i}
//
// `local_type` is the per-site number operator: "N" for the spinless fermion /
// boson blocks, "Ntot" for the spinful electron block.
OpRule totaln_expansion_rule(int64_t nsites, std::string const &local_type = "N");

} // namespace xdiag::algebra
