#include "algebra.hpp"

#include <xdiag/algebra/apply.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/parallel/mpi/cdot_distributed.hpp>
#endif

namespace xdiag {

double norm(State const &v) try {
  if (v.n_cols() > 1) {
    XDIAG_THROW("Cannot compute norm of state with more than one column");
    return 0;
  } else {
    if (v.isreal()) {
      return norm(v.block(), v.vector(0, false));
    } else {
      return norm(v.block(), v.vectorC(0, false));
    }
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double dot(State const &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
    return 0;
  }

  if ((v.n_cols() > 1) || (w.n_cols() > 1)) {
    XDIAG_THROW(
        "Cannot compute dot product of state with more than one column");
  }

  if ((v.isreal()) && (w.isreal())) {
    return dot(v.block(), v.vector(0, false), w.vector(0, false));
  } else {
    XDIAG_THROW("Unable to compute real dot product of a complex "
                "state, consider using dotC instead.");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

complex dotC(State const &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
  }

  if ((v.n_cols() > 1) || (w.n_cols() > 1)) {
    XDIAG_THROW(
        "Cannot compute dotC product of state with more than one column");
  }
  if ((v.isreal()) && (w.isreal())) {
    return dot(v.block(), v.vector(0, false), w.vector(0, false));
  } else if ((v.isreal()) && (!w.isreal())) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v2.vectorC(0, false), w.vectorC(0, false));
  } else if ((w.isreal()) && (!v.isreal())) {
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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

double inner(OpSum const &ops, State const &v) try {
  if (v.isreal() && ops.isreal()) {
    auto w = v;
    apply(ops, v, w);
    return dot(w, v);
  } else {
    XDIAG_THROW("\"inner\" function computing product <psi | O | psi> can only "
                "be called if both the state and the Ops are real. Maybe use "
                "innerC(...) instead.");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

double inner(Op const &op, State const &v) try {
  return inner(OpSum({op}), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

complex innerC(OpSum const &ops, State const &v) try {
  if (v.isreal() && ops.isreal()) {
    auto w = v;
    apply(ops, v, w);
    return (complex)dot(w, v);
  } else if (v.isreal() && !ops.isreal()) {
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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

complex innerC(Op const &op, State const &v) try {
  return innerC(OpSum({op}), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v) try {
  return std::sqrt(xdiag::real(dot(block, v, v)));
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template double norminf(Block const &, arma::Col<double> const &);
template double norminf(Block const &, arma::Col<complex> const &);

State &operator*=(State &X, double alpha) try {
  if (X.isreal()) {
    X.matrix(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(State const &X, double alpha) try {
  auto res = X;
  res *= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(double alpha, State const &X) try {
  return X * alpha;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator*=(State &X, complex alpha) try {
  if (X.isreal()) {
    X.make_complex();
    X.matrixC(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(State const &X, complex alpha) try {
  auto res = X;
  res *= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(complex alpha, State const &X) try {
  return X * alpha;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator/=(State &X, double alpha) try {
  if (X.isreal()) {
    X.matrix(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator/(State const &X, double alpha) try {
  auto res = X;
  res /= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator/=(State &X, complex alpha) try {
  if (X.isreal()) {
    X.make_complex();
    X.matrixC(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator/(State const &X, complex alpha) try {
  auto res = X;
  res /= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator+=(State &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot add two states from different blocks");
  }

  if (v.isreal() && w.isreal()) {
    v.matrix(false) += w.matrix(false);
  } else if (!v.isreal() && w.isreal()) {
    auto wcplx = w;
    wcplx.make_complex();
    v.matrixC(false) += wcplx.matrixC(false);
  } else if (v.isreal() && !w.isreal()) {
    v.make_complex();
    v.matrixC(false) += w.matrixC(false);
  } else { // (!v.isreal() && !w.isreal())
    v.matrixC(false) += w.matrixC(false);
  }

  return v;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator-=(State &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot subtract two states from different blocks");
  }

  if (v.isreal() && w.isreal()) {
    v.matrix(false) -= w.matrix(false);
  } else if (!v.isreal() && w.isreal()) {
    auto wcplx = w;
    wcplx.make_complex();
    v.matrixC(false) -= wcplx.matrixC(false);
  } else if (v.isreal() && !w.isreal()) {
    v.make_complex();
    v.matrixC(false) -= w.matrixC(false);
  } else { // (!v.isreal() && !w.isreal())
    v.matrixC(false) -= w.matrixC(false);
  }

  return v;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator+(State const &v, State const &w) try {
  auto res = v;
  res += w;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator-(State const &v, State const &w) try {
  auto res = v;
  res -= w;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag
