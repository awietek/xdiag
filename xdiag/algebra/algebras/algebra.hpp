// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/rules.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>

namespace xdiag::algebra {

// Algebra describes how to bring an OpSum into normal order for a specific
// physical system. It contains:
//
//   fermionic_types   — the Op types that anticommute at different sites
//                       (used to determine swap signs)
//   allowed_types     — all valid input Op types (checked before expansion)
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
  int64_t nsites;   // number of sites 
  int64_t d;        // local Hilbert space dimension per site
  std::set<std::string> fermionic_types;
  std::set<std::string>
      allowed_types; // all valid input types (checked before expansion)
  std::vector<OpRule> expansion_rules;
  std::vector<MonomialRule> algebra_rules;
};

Algebra implementation_algebra(Boson const &block);
Algebra implementation_algebra(Spinhalf const &block);
Algebra implementation_algebra(Fermion const &block);
Algebra implementation_algebra(Electron const &block);
Algebra implementation_algebra(tJ const &block);
Algebra implementation_algebra(Block const &block);

Algebra symmetry_algebra(Boson const &block);
Algebra symmetry_algebra(Spinhalf const &block);
Algebra symmetry_algebra(Fermion const &block);
Algebra symmetry_algebra(Electron const &block);
Algebra symmetry_algebra(tJ const &block);
Algebra symmetry_algebra(Block const &block);

} // namespace xdiag::algebra
