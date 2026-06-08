// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <variant>

#include <xdiag/algebra/algebras/matrix_algebra.hpp>
#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::algebra {

Algebra implementation_algebra(Spinhalf const &) try {
  return spinhalf_implementation_algebra();
}
XDIAG_CATCH

Algebra implementation_algebra(Boson const &block) try {
  return matrix_algebra(block.d());
}
XDIAG_CATCH

Algebra implementation_algebra(Block const &block) try {
  return std::visit([](auto const &b) { return implementation_algebra(b); },
                    block);
}
XDIAG_CATCH

Algebra symmetry_algebra(Spinhalf const &) try { return matrix_algebra(2); }
XDIAG_CATCH

Algebra symmetry_algebra(Boson const &block) try {
  return matrix_algebra(block.d());
}
XDIAG_CATCH

Algebra symmetry_algebra(Block const &block) try {
  return std::visit([](auto const &b) { return symmetry_algebra(b); }, block);
}
XDIAG_CATCH

} // namespace xdiag::algebra
