// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

struct Algebra;

// Generic site-major sorting step shared by the symmetry algebras: swap two
// adjacent single-site operators on different sites into ascending Op order,
// with a (-1) sign when both are fermionic. nullopt if already sorted.
std::optional<OpSum> sort_sites_step(Monomial const &mono, Algebra const &alg);

} // namespace xdiag::algebra
