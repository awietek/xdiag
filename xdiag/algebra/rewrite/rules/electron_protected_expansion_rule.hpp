// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <set>
#include <string>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule that expands the named composite electron operators (Hopup,
// Hopdn, Hop*Asym, Nup, Ndn, Nupdn, N*N*, SzSz, HubbardU) into the elementary
// generators Cdagup/Cup/Cdagdn/Cdn. The spinful analogue of
// fermion_protected_expansion_rule.
//
// A type listed in `protected_types` is left untouched when it is the sole
// operator of a monomial (size 1), so named operators with a dedicated matrix
// kernel survive and dispatch directly. As soon as such an operator appears
// inside a product (size > 1) it is expanded, so the whole product can be
// normal ordered into a Cdag/C string (handled by term_cdagc_string).
//
// Repeated application to a fixed point expands every composite operator that
// occurs in a product. Operators NOT in `protected_types` are always expanded.
MonomialRule
electron_protected_expansion_rule(std::set<std::string> const &protected_types,
                                  int64_t nsites);

} // namespace xdiag::algebra
