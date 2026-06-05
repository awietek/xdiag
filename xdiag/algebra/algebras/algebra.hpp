// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/rules.hpp>
#include <xdiag/blocks/blocks.hpp>

namespace xdiag::algebra {

// Algebra describes how to bring an OpSum into normal order for a specific
// physical system. It contains:
//
//   elementary_types  — the canonical generators (e.g. {S+, S-, Sz} for spin)
//   fermionic_types   — subset of elementary_types that anticommute at
//                       different sites (used to determine swap signs)
//   expansion_rules   — OpRules that expand compound operators into products
//                       of elementary ones (applied first, to fixed point)
//   algebra_rules     — MonomialRules that simplify same-site products and
//                       sort operators into canonical order (applied after
//                       expansion, to fixed point)
//
// Concrete algebras are defined in their own headers (spin_algebra.hpp,
// electron_algebra.hpp, tj_algebra.hpp, matrix_algebra.hpp,
// spinhalf_implementation_algebra.hpp).
struct Algebra {
  std::string name; // human-readable name for errors
  int64_t d;        // local Hilbert space dimension per site
  std::vector<std::string> elementary_types;
  std::set<std::string> fermionic_types;
  std::set<std::string>
      allowed_types; // all valid input types (checked before expansion)
  std::vector<OpRule> expansion_rules;
  std::vector<MonomialRule> algebra_rules;
};

// Returns the Algebra used to bring an OpSum into normal order for the given
// block. Each block type has a canonical algebra (e.g. Spinhalf uses the
// spin-1/2 implementation algebra). Used wherever an OpSum must be normal
// ordered or compared (e.g. isapprox / hermiticity checks) in the context of a
// concrete block.
Algebra algebra(Spinhalf const &block);
Algebra algebra(Block const &block);

} // namespace xdiag::algebra
