// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule: spinless-fermion same-site CAR (adjacent pair, same site):
//   Cdag*Cdag = 0
//   C*C       = 0
//   C*Cdag    = Id - Cdag*C        ({C, Cdag} = 1)
//   Cdag*C    = (already canonical, left untouched)
// Canonical same-site order places the creation operator on the left
// (Cdag*C). Only the out-of-order / nilpotent directions fire.
MonomialRule fermion_same_site_rule();

} // namespace xdiag::algebra
