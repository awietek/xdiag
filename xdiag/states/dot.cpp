// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "dot.hpp"

#include <xdiag/math/dot.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {
double dot(State const &v, State const &w) try {
  if ((!isvalid(v)) || (!isvalid(w))) {
    return 0.;
  }

  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
    return 0;
  }

  if ((v.ncols() > 1) || (w.ncols() > 1)) {
    XDIAG_THROW("Cannot compute dot product of state with more than one "
                "column. Consider using matrix_dot instead.");
  }

  if ((isreal(v)) && (isreal(w))) {
    return math::dot(v.block(), v.vector(0, false), w.vector(0, false));
  } else {
    XDIAG_THROW("Unable to compute real dot product of a complex "
                "state. Consider using dotC instead.");
  }
}
XDIAG_CATCH

complex dotC(State const &v, State const &w) try {
  if ((!isvalid(v)) || (!isvalid(w))) {
    return complex(0.);
  }

  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
  }

  if ((v.ncols() > 1) || (w.ncols() > 1)) {
    XDIAG_THROW("Cannot compute dotC product of state with more than one "
                "column. Consider using matrix_dotC instead.");
  }
  if ((isreal(v)) && (isreal(w))) {
    return math::dot(v.block(), v.vector(0, false), w.vector(0, false));
  } else if ((isreal(v)) && (!isreal(w))) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return math::dot(v.block(), v2.vectorC(0, false), w.vectorC(0, false));
  } else if ((isreal(w)) && (!isreal(v))) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return math::dot(v.block(), v.vectorC(0, false), w2.vectorC(0, false));
  } else {
    return math::dot(v.block(), v.vectorC(0, false), w.vectorC(0, false));
  }
}
XDIAG_CATCH

arma::mat matrix_dot(State const &v, State const &w) try {
  if ((!isvalid(v)) || (!isvalid(w))) {
    return arma::mat();
  }

  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
  }

  if ((isreal(v)) && (isreal(w))) {
    return math::matrix_dot(v.block(), v.matrix(false), w.matrix(false));
  } else {
    XDIAG_THROW("Unable to compute real dot product of a complex "
                "state. Consider using dotC instead.");
  }
}
XDIAG_CATCH

arma::cx_mat matrix_dotC(State const &v, State const &w) try {
  if ((!isvalid(v)) || (!isvalid(w))) {
    return arma::cx_mat();
  }

  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
  }

  if ((isreal(v)) && (isreal(w))) {
    return arma::conv_to<arma::cx_mat>::from(
        math::matrix_dot(v.block(), v.matrix(false), w.matrix(false)));
  } else if ((isreal(v)) && (!isreal(w))) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return math::matrix_dot(v.block(), v2.matrixC(false), w.matrixC(false));
  } else if ((isreal(w)) && (!isreal(v))) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return math::matrix_dot(v.block(), v.matrixC(false), w2.matrixC(false));
  } else {
    return math::matrix_dot(v.block(), v.matrixC(false), w.matrixC(false));
  }
}
XDIAG_CATCH

} // namespace xdiag
