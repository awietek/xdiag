// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <xdiag/algebra/apply.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/utils/arma_to_cx.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/parallel/mpi/cdot_distributed.hpp>
#endif

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
      return norm(v.block(), v.vector(0, false));
    } else {
      return norm(v.block(), v.vectorC(0, false));
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
      return norm1(v.block(), v.vector(0, false));
    } else {
      return norm1(v.block(), v.vectorC(0, false));
    }
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double norminf(State const &v) try {
  if (!isvalid(v)) {
    return 0.;
  }

  if (v.ncols() > 1) {
    XDIAG_THROW("Cannot compute norm of state with more than one column");
    return 0;
  } else {
    if (isreal(v)) {
      return norminf(v.block(), v.vector(0, false));
    } else {
      return norminf(v.block(), v.vectorC(0, false));
    }
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

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
    return dot(v.block(), v.vector(0, false), w.vector(0, false));
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
    return dot(v.block(), v.vector(0, false), w.vector(0, false));
  } else if ((isreal(v)) && (!isreal(w))) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v2.vectorC(0, false), w.vectorC(0, false));
  } else if ((isreal(w)) && (!isreal(v))) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v.vectorC(0, false), w2.vectorC(0, false));
  } else {
    return dot(v.block(), v.vectorC(0, false), w.vectorC(0, false));
  }
}
XDIAG_CATCH

arma::mat matrix_dot(State const &v, State const &w) try {
  if ((!isvalid(v)) || (!isvalid(w))) {
    return arma::mat();
  }

  if (v.block() != w.block()) {
    XDIAG_THROW(
        "Cannot form matrix_dot product for states on different blocks");
    return 0;
  }

  if ((isreal(v)) && (isreal(w))) {
    return matrix_dot(v.block(), v.matrix(false), w.matrix(false));
  } else {
    XDIAG_THROW("Unable to compute real matrix_dot product of a complex "
                "state. Consider using matrix_dotC instead.");
  }
}
XDIAG_CATCH

arma::cx_mat matrix_dotC(State const &v, State const &w) try {
  if ((!isvalid(v)) || (!isvalid(w))) {
    return arma::cx_mat();
  }

  if (v.block() != w.block()) {
    XDIAG_THROW(
        "Cannot form matrix_dot product for states on different blocks");
    return 0;
  }

  if (isreal(v) && isreal(w)) {
    return utils::to_cx_mat(
        matrix_dot(v.block(), v.matrix(false), w.matrix(false)));
  } else if ((isreal(v)) && (!isreal(w))) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return matrix_dot(v.block(), v2.matrixC(false), w.matrixC(false));
  } else if ((isreal(w)) && (!isreal(v))) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return matrix_dot(v.block(), v.matrixC(false), w2.matrixC(false));
  } else {
    return matrix_dot(v.block(), v.matrixC(false), w.matrixC(false));
  }
}
XDIAG_CATCH

double inner(OpSum const &ops, State const &v) try {
  if (!isvalid(v)) {
    return 0.;
  }

  if (isreal(v) && isreal(ops)) {
    auto w = v;
    apply(ops, v, w);
    return dot(w, v);
  } else {
    XDIAG_THROW("\"inner\" function computing product <psi | O | psi> can only "
                "be called if both the state and the Ops are real. Maybe use "
                "innerC(...) instead.");
  }
}
XDIAG_CATCH

double inner(Op const &op, State const &v) try {
  if (!isvalid(v)) {
    return 0.;
  }

  return inner(OpSum(op), v);
}
XDIAG_CATCH

complex innerC(OpSum const &ops, State const &v) try {
  if (!isvalid(v)) {
    return complex(0.);
  }

  if (isreal(v) && isreal(ops)) {
    auto w = v;
    apply(ops, v, w);
    return (complex)dot(w, v);
  } else if (isreal(v) && !isreal(ops)) {
    auto v2 = v;
    auto w = v2;
    v2.make_complex();
    apply(ops, v2, w);
    return dotC(w, v);
  } else {
    auto w = v;
    apply(ops, v, w);
    return dotC(w, v);
  }
}
XDIAG_CATCH

complex innerC(Op const &op, State const &v) try {
  if (!isvalid(v)) {
    return complex(0.);
  }
  return innerC(OpSum(op), v);
}
XDIAG_CATCH

double dot(Block const &block, arma::vec const &v, arma::vec const &w) try {
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return arma::dot(v, w);
#ifdef XDIAG_USE_MPI
  }
#endif
}
XDIAG_CATCH

complex dot(Block const &block, arma::cx_vec const &v,
            arma::cx_vec const &w) try {
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return arma::cdot(v, w);
#ifdef XDIAG_USE_MPI
  }
#endif
}
XDIAG_CATCH

template <typename coeff_t>
arma::Mat<coeff_t> matrix_dot(Block const &block, arma::Mat<coeff_t> const &V,
                              arma::Mat<coeff_t> const &W) try {
  if (V.n_rows != V.n_rows) {
    XDIAG_THROW("Input matrices do not have the same number of rows");
  }
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    int64_t L = V.n_rows;
    int64_t m = V.n_cols;
    int64_t n = W.n_cols;
    arma::Mat<coeff_t> result(m, n, arma::fill::zeros);
    for (int64_t i = 0; i < m; ++i) {
      for (int64_t j = 0; j < n; ++j) {
        arma::Col<coeff_t> cv(const_cast<coeff_t *>(V.colptr(i)), L, false,
                              true);
        arma::Col<coeff_t> cw(const_cast<coeff_t *>(W.colptr(j)), L, false,
                              true);
        result(i, j) = cdot_distributed(cv, cw);
      }
    }
    return result;
  } else {
#else
  (void)block;
#endif
    return V.t() * W;
#ifdef XDIAG_USE_MPI
  }
#endif
}
XDIAG_CATCH

template arma::mat matrix_dot(Block const &, arma::mat const &,
                              arma::mat const &);
template arma::cx_mat matrix_dot(Block const &, arma::cx_mat const &,
                                 arma::cx_mat const &);

template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v) try {
  return std::sqrt(xdiag::real(dot(block, v, v)));
}
XDIAG_CATCH

template double norm(Block const &, arma::Col<double> const &);
template double norm(Block const &, arma::Col<complex> const &);

template <typename coeff_t>
double norm1(Block const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, 1);
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return nrm;
}
XDIAG_CATCH

template double norm1(Block const &, arma::Col<double> const &);
template double norm1(Block const &, arma::Col<complex> const &);

template <typename coeff_t>
double norminf(Block const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, "inf");
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  return nrm;
}
XDIAG_CATCH

template double norminf(Block const &, arma::Col<double> const &);
template double norminf(Block const &, arma::Col<complex> const &);

} // namespace xdiag
