// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <variant>

#include <xdiag/algebra/algebras/fermion_algebra.hpp>
#include <xdiag/algebra/algebras/fermion_implementation_algebra.hpp>
#include <xdiag/algebra/algebras/matrix_algebra.hpp>
#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag::algebra {

Algebra implementation_algebra(Spinhalf const &block) try {
  return spinhalf_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(Boson const &block) try {
  return matrix_algebra(block.nsites(), block.d());
}
XDIAG_CATCH

Algebra implementation_algebra(Fermion const &block) try {
  return fermion_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(Block const &block) try {
  // Exhaustive per-type dispatch: adding a new alternative to the Block variant
  // without a matching overload here is a compile-time error (std::visit can no
  // longer call the visitor for every alternative).
  return std::visit(
      utils::overload{
          [](Spinhalf const &b) { return implementation_algebra(b); },
          [](Boson const &b) { return implementation_algebra(b); },
          [](Fermion const &b) { return implementation_algebra(b); },
      },
      block);
}
XDIAG_CATCH

Algebra symmetry_algebra(Spinhalf const &block) try {
  return matrix_algebra(block.nsites(), 2);
}
XDIAG_CATCH

Algebra symmetry_algebra(Boson const &block) try {
  return matrix_algebra(block.nsites(), block.d());
}
XDIAG_CATCH

Algebra symmetry_algebra(Fermion const &block) try {
  return fermion_algebra(block.nsites());
}
XDIAG_CATCH

Algebra symmetry_algebra(Block const &block) try {
  // Exhaustive per-type dispatch (see implementation_algebra(Block) above).
  return std::visit(
      utils::overload{
          [](Spinhalf const &b) { return symmetry_algebra(b); },
          [](Boson const &b) { return symmetry_algebra(b); },
          [](Fermion const &b) { return symmetry_algebra(b); },
      },
      block);
}
XDIAG_CATCH

} // namespace xdiag::algebra
