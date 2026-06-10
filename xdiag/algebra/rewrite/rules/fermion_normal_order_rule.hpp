// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule: brings a product of elementary spinless-fermion operators into
// textbook normal order:
//
//   all Cdag (ascending site) to the left of all C (ascending site)
//
// i.e. canonical order is: creation before annihilation, then ascending site.
// Every adjacent swap of two (different-site) fermionic operators contributes a
// -1 sign. Same-site pairs are left to fermion_same_site_rule (this rule skips
// them). Non Cdag/C operators are ignored.
//
// Unlike sort_sites_rule (which orders purely by site and interleaves Cdag/C),
// this produces the creation-left form expected by term_cdagc_string.
MonomialRule fermion_normal_order_rule();

} // namespace xdiag::algebra
