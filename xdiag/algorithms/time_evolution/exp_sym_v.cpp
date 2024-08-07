#include "exp_sym_v.hpp"

#include <xdiag/algebra/apply.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>

namespace xdiag {

template <typename coeff_t, class multiply_f>
double exp_sym_v_inplace(multiply_f mult, arma::Col<coeff_t> &X, coeff_t tau,
                         bool normalize = false, double shift = 0,
                         double precision = 1e-12,
                         int64_t max_iterations = 1000,
                         double deflation_tol = 1e-7) try {

  double norm = arma::norm(X);
  auto v0 = X;

  auto dot = [](arma::Col<coeff_t> const &v, arma::Col<coeff_t> const &w) {
    return arma::cdot(v, w);
  };
  auto converged = [precision, tau, norm](Tmatrix const &tmat) {
    return lanczos::converged_time_evolution(tmat, tau, precision, norm);
  };

  auto operation_void = [](arma::Col<coeff_t> const &) {};
  auto r = lanczos::lanczos(mult, dot, converged, operation_void, v0,
                            max_iterations, deflation_tol);

  double e0 = r.eigenvalues(0);

  // Compute the tridiagonal matrix
  arma::mat tmat = arma::diagmat(r.alphas);
  if (r.alphas.n_rows > 1) {
    tmat += arma::diagmat(r.betas.head(r.betas.size() - 1), 1) +
            arma::diagmat(r.betas.head(r.betas.size() - 1), -1);
  }

  // Subtract shift from diagonal
  if (shift != 0.) {
    for (int64_t i = 0; i < (int64_t)tmat.n_cols; ++i) {
      tmat(i, i) -= shift;
    }
  }

  arma::Mat<coeff_t> texp = arma::expmat(tau * tmat);
  arma::Col<coeff_t> linear_combination = texp.col(0);

  v0 = X;
  X.zeros();
  int64_t iter = 0;
  auto mult2 = [&](arma::Col<coeff_t> const &v, arma::Col<coeff_t> &w) {
    mult(v, w);
    ++iter;
  };
  auto operation = [&linear_combination, &iter,
                    &X](arma::Col<coeff_t> const &v) {
    X += linear_combination(iter) * v;
  };

  lanczos::lanczos(mult2, dot, converged, operation, v0, max_iterations,
                   deflation_tol);

  if (!normalize) {
    X *= norm;
  } else {
    double nrm = arma::norm(X);
    X /= nrm;
  }

  return e0;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return 0.;
}

State exp_sym_v(OpSum const &ops, State state, double tau, bool normalize,
                double shift, double precision, int64_t max_iterations,
                double deflation_tol) try {
  auto const &block = state.block();

  // Real time evolution is possible
  if (state.isreal() && ops.isreal()) {
    int iter = 1;
    auto mult = [&iter, &ops, &block](arma::vec const &v, arma::vec &w) {
      auto ta = rightnow();
      apply(ops, block, v, block, w);
      Log(2, "Lanczos iteration {}", iter);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };
    arma::vec v = state.vector(0, false);
    exp_sym_v_inplace(mult, v, tau, normalize, shift, precision, max_iterations,
                      deflation_tol);
    return state;

    // Refer to complex time evolution
  } else {
    return exp_sym_v(ops, state, complex(tau), normalize, shift, precision,
                     max_iterations, deflation_tol);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return State();
}

State exp_sym_v(OpSum const &ops, State state, complex tau, bool normalize,
                double shift, double precision, int64_t max_iterations,
                double deflation_tol) try {
  auto const &block = state.block();
  state.make_complex();

  int iter = 1;
  auto mult = [&iter, &ops, &block](arma::cx_vec const &v, arma::cx_vec &w) {
    auto ta = rightnow();
    apply(ops, block, v, block, w);
    Log(2, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };
  arma::cx_vec v = state.vectorC(0, false);
  exp_sym_v_inplace(mult, v, tau, normalize, shift, precision, max_iterations,
                    deflation_tol);
  return state;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return State();
}
} // namespace xdiag
