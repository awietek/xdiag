// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <set>
#include <string>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule that expands the named composite tJ operators (Hopup, Hopdn,
// Hop*Asym, Nup, Ndn, N*N*, SzSz, tJSzSz, Exchange) into the elementary
// generators Cdagup/Cup/Cdagdn/Cdn. The tJ analogue of
// electron_protected_expansion_rule -- it additionally knows the tJ-specific
// Exchange and tJSzSz composites.
//
// A type listed in `protected_types` is left untouched when it is the sole
// operator of a monomial (size 1), so named operators with a dedicated matrix
// kernel survive and dispatch directly. As soon as such an operator appears
// inside a product (size > 1) it is expanded, so the whole product can be
// normal ordered into a Cdag/C string (handled by term_cdagc_string).
MonomialRule
tj_protected_expansion_rule(std::set<std::string> const &protected_types,
                            int64_t nsites);

} // namespace xdiag::algebra
