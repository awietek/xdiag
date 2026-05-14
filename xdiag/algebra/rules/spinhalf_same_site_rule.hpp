// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// MonomialRule: spin-1/2 same-site relations (adjacent pair, same site):
//   S+*S+ = 0         S-*S- = 0
//   S+*Sz = -1/2 S+   Sz*S+ = +1/2 S+
//   S-*Sz = +1/2 S-   Sz*S- = -1/2 S-
//   S+*S- = 1/2 Id + Sz
//   S-*S+ = 1/2 Id - Sz
//   Sz*Sz = 1/4 Id
MonomialRule spinhalf_same_site_rule();

} // namespace xdiag::operators
