#pragma once

#include <hydra/algorithms/lanczos/lanczos.h>
#include <hydra/algorithms/lanczos/lanczos_build.h>
#include <hydra/algorithms/lanczos/lanczos_convergence.h>
#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/utils/timing.h>

namespace hydra {

template <typename coeff_t, class multiply_f>
double exp_sym_v_inplace(multiply_f A, arma::Col<coeff_t> &X, coeff_t tau,
                         double precision = 1e-12, bool shift = false) {

  double norm = arma::norm(X);

  auto converged = [precision, tau, norm](Tmatrix const &tmat) {
    return converged_time_evolution(tmat, tau, precision, norm);
  };
  auto v0 = X;
  auto tmat = lanczos(A, v0, converged);
  double e0 = tmat.eigenvalues()(0);

  // Cast to complex (TODO: make this generic)
  arma::mat tmatc = tmat.mat();

  if (shift) {
    for (arma::uword i = 0; i < tmat.size(); ++i) {
      tmatc(i, i) -= e0;
    }
  }

  arma::Mat<coeff_t> texp = arma::expmat(tau * tmatc);
  arma::Col<coeff_t> linear_combination = texp.col(0);

  auto [tmat2, vec] = lanczos_build(A, X, linear_combination);
  X = vec * norm;
  return e0;
}

template <typename coeff_t, class multiply_f>
arma::Col<coeff_t> exp_sym_v(multiply_f A, arma::Col<coeff_t> const &X,
                             coeff_t tau, double precision = 1e-12) {
  auto v0 = X;
  exp_sym_v_inplace(A, v0, tau, precision);
  return v0;
}

template <class StateT>
StateT exp_sym_v(BondList const &bonds, StateT const &state,
                 typename StateT::coeff_t tau, double precision = 1e-12) {
  auto const &block = state.block();

  using coeff_t = typename StateT::coeff_t;

  int iter = 1;
  auto mult = [&iter, &bonds, &block](arma::Col<coeff_t> const &v,
                                      arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    apply(bonds, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };

  auto v0 = state.vector();
  auto t0 = rightnow();
  exp_sym_v_inplace(mult, v0, tau, precision);
  timing(t0, rightnow(), "Lanczos time", 1);
  return StateT(block, v0);
}

} // namespace hydra
