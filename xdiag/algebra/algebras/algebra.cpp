// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <variant>

#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::algebra {

Algebra algebra(Spinhalf const &) try {
  return spinhalf_implementation_algebra();
}
XDIAG_CATCH

Algebra algebra(Block const &block) try {
  return std::visit([](auto const &b) { return algebra(b); }, block);
}
XDIAG_CATCH

} // namespace xdiag::algebra
