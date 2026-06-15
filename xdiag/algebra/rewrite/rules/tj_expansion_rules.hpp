// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>
#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// tJ expansion rules: all electron_expansion_rules plus the tJ-specific
// tJSzSz, tJSdotS, and (in tJ) Exchange expansions in terms of c/c†.
std::vector<OpRule> tj_expansion_rules();

} // namespace xdiag::algebra
