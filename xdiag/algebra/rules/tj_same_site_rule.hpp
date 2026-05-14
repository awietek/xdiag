// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// MonomialRule: tJ same-site CAR (projected, no double occupancy):
//   Cdagup*Cdagup = 0    Cup*Cup = 0
//   Cdagdn*Cdagdn = 0    Cdn*Cdn = 0
//   Cup*Cdagup    = Id - Cdagup*Cup - Cdagdn*Cdn   ({Cdagup,Cup} = Id - Ndn)
//   Cdn*Cdagdn    = Id - Cdagdn*Cdn - Cdagup*Cup   ({Cdagdn,Cdn} = Id - Nup)
//   Cdn*Cdagup    = 0   (no double occupancy)
//   Cup*Cdagdn    = 0   (no double occupancy)
//   Cdagup*Cdagdn = -Cdagdn*Cdagup
//   Cup*Cdn       = -Cdn*Cup
//   (only out-of-canonical-order direction fires; canonical order is
//    alphabetical Cdagdn < Cdagup < Cdn < Cup)
MonomialRule tj_same_site_rule();

} // namespace xdiag::operators
