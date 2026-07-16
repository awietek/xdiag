// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "norm.hpp"

#include <xdiag/math/norm.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

double norm(State const &v) try {
  if (!isvalid(v)) {
    return 0.;
  }

  if (v.ncols() > 1) {
    XDIAG_THROW("Cannot compute norm of state with more than one column");
    return 0;
  } else {
    if (isreal(v)) {
      return math::norm(v.block(), v.vector(0, false));
    } else {
      return math::norm(v.block(), v.vectorC(0, false));
    }
  }
}
XDIAG_CATCH

double norm1(State const &v) try {
  if (!isvalid(v)) {
    return 0.;
  }

  if (v.ncols() > 1) {
    XDIAG_THROW("Cannot compute norm of state with more than one column");
    return 0;
  } else {
    if (isreal(v)) {
      return math::norm1(v.block(), v.vector(0, false));
    } else {
      return math::norm1(v.block(), v.vectorC(0, false));
    }
  }
}
XDIAG_CATCH

double norminf(State const &v) try {
  if (!isvalid(v)) {
    return 0.;
  }

  if (v.ncols() > 1) {
    XDIAG_THROW("Cannot compute norm of state with more than one column");
    return 0;
  } else {
    if (isreal(v)) {
      return math::norminf(v.block(), v.vector(0, false));
    } else {
      return math::norminf(v.block(), v.vectorC(0, false));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag
