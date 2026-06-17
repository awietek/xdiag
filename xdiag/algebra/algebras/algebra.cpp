// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <variant>

#include <xdiag/algebra/algebras/electron_algebra.hpp>
#include <xdiag/algebra/algebras/electron_implementation_algebra.hpp>
#include <xdiag/algebra/algebras/fermion_algebra.hpp>
#include <xdiag/algebra/algebras/fermion_implementation_algebra.hpp>
#include <xdiag/algebra/algebras/matrix_algebra.hpp>
#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/algebra/algebras/tj_algebra.hpp>
#include <xdiag/algebra/algebras/tj_implementation_algebra.hpp>
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

Algebra implementation_algebra(Electron const &block) try {
  return electron_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(tJ const &block) try {
  return tj_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(Block const &block) try {
  // Generic visitor: a block type without an implementation_algebra overload is
  // a compile error when the lambda is instantiated for that alternative.
  return std::visit(
      [](auto const &b) { return implementation_algebra(b); }, block);
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

Algebra symmetry_algebra(Electron const &block) try {
  return electron_algebra(block.nsites());
}
XDIAG_CATCH

Algebra symmetry_algebra(tJ const &block) try {
  return tj_algebra(block.nsites());
}
XDIAG_CATCH

Algebra symmetry_algebra(Block const &block) try {
  // Generic visitor: a block type without a symmetry_algebra overload is a
  // compile error when the lambda is instantiated for that alternative.
  return std::visit([](auto const &b) { return symmetry_algebra(b); }, block);
}
XDIAG_CATCH

} // namespace xdiag::algebra
