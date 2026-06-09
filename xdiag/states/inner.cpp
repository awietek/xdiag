// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "inner.hpp"

#include <xdiag/states/apply.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

double inner(State const &v, Op const &A, State const &w) try {
  return inner(v, OpSum(A), w);
}
XDIAG_CATCH
double inner(State const &v, Monomial const &A, State const &w) try {
  return inner(v, OpSum(A), w);
}
XDIAG_CATCH
double inner(State const &v, OpSum const &A, State const &w) try {
  if (isreal(v) && isreal(A) && isreal(w)) {
    return dot(v, apply(A, w));
  } else {
    XDIAG_THROW(
        "Both the input states and the OpSum need to be real for inner(...). "
        "Consider using innerC(...) for complex arguments instead");
  }
}
XDIAG_CATCH

complex innerC(State const &v, Op const &A, State const &w) try {
  return innerC(v, OpSum(A), w);
}
XDIAG_CATCH
complex innerC(State const &v, Monomial const &A, State const &w) try {
  return innerC(v, OpSum(A), w);
}
XDIAG_CATCH
complex innerC(State const &v, OpSum const &A, State const &w) try {
  return dotC(v, apply(A, w));
}
XDIAG_CATCH

} // namespace xdiag
