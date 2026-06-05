// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule: electron same-site CAR (adjacent pair, same site):
//   Cdagup*Cdagup = 0    Cup*Cup = 0
//   Cdagdn*Cdagdn = 0    Cdn*Cdn = 0
//   Cup*Cdagup    = Id - Cdagup*Cup
//   Cdn*Cdagdn    = Id - Cdagdn*Cdn
//   Cdagup*Cdagdn = -Cdagdn*Cdagup
//   Cup*Cdn       = -Cdn*Cup
//   Cdagup*Cdn    = -Cdn*Cdagup
//   Cup*Cdagdn    = -Cdagdn*Cup
//   (only out-of-canonical-order direction fires; canonical order is
//    alphabetical Cdagdn < Cdagup < Cdn < Cup)
MonomialRule electron_same_site_rule();

} // namespace xdiag::algebra
