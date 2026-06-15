// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// Normal-ordering rule for strings of elementary electron operators
// (Cdagup/Cup/Cdagdn/Cdn) used by the electron IMPLEMENTATION algebra. It brings
// a product into the canonical order expected by the matrix kernel
// (matrices/blocks/electron/terms/term_cdagc_string.hpp), namely the
// "all ups, then all dns" Jordan-Wigner order of the electron basis:
//
//   * all up operators (Cdagup/Cup) to the LEFT of all dn operators
//     (Cdagdn/Cdn) -- different fermionic modes anticommute, sign -1 per swap;
//   * within each sector, creation before annihilation, then ascending site;
//   * same (sector, site) pairs are resolved by the canonical anticommutation
//     relations: Cdag*Cdag = C*C = 0, Cdag*C is kept (number operator),
//     C*Cdag = Id - Cdag*C.
//
// NOTE: this differs from electron_same_site_rule / the symmetry algebra, which
// use a creation-left, dn-before-up convention suited to quantum-number
// analysis. Here the convention must match the basis Jordan-Wigner ordering.
MonomialRule electron_normal_order_rule();

} // namespace xdiag::algebra
