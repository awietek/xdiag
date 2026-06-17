// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule that brings a product of elementary tJ operators
// (Cdagup/Cup/Cdagdn/Cdn) into the "all ups then all dns, creation-left,
// ascending" normal order expected by term_cdagc_string.
//
// Different sites anticommute like ordinary fermions (swap with a -1 sign).
// Same-site pairs use the tJ (projected) relations rather than the canonical
// ones -- the key differences from the free electron normal ordering are:
//   Cup_i Cdagup_i  = Id - Cdagup_i Cup_i - Cdagdn_i Cdn_i   ({Cup,Cdagup}=1-Ndn)
//   Cdn_i Cdagdn_i  = Id - Cdagdn_i Cdn_i - Cdagup_i Cup_i
//   Cup_i Cdagdn_i  = 0,   Cdn_i Cdagup_i = 0                (no double occupancy)
//   Cdagdn_i Cup_i  (S-) is kept as-is (NOT swapped to Cup_i Cdagdn_i, which
//                    would be zero); Cdagup_i Cdn_i (S+) likewise.
// Double creations / annihilations at a site (e.g. Cdagup_i Cdagdn_i) are merely
// reordered here; term_cdagc_string's no-double-occupancy projection sends them
// to zero when applied.
MonomialRule tj_normal_order_rule();

} // namespace xdiag::algebra
