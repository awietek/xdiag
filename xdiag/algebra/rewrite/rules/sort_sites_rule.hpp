// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>
#include <string>
#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule: sort adjacent single-site operators by site index.
// fermionic_types names the Op types that anticommute when on different
// sites; swapping a pair of fermionic ops contributes a -1 sign.
MonomialRule sort_sites_rule(std::set<std::string> const &fermionic_types);

} // namespace xdiag::algebra
