//
// Created by Alex Wietek on 7/3/23
//
#pragma once

#include <string>

#include "extern/armadillo/armadillo"
#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra {

template <typename coeff_t> struct lanczos_result {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  arma::Mat<coeff_t> V;
  arma::Col<coeff_t> residual;
  int num_iterations;
  std::string criterion;
  double anorm;
  arma::vec omegas;
}

// Lanczos implementation with partial reorthogonalization

// implementation guided by
// https://github.com/areslp/matlab/blob/master/PROPACK/lanpro.m
template <class coeff_t, class multiply_f, class convergence_f>
lanczos_result<coeff_t>
lanczos_pro(multiply_f mult, arma::Col<coeff_t> &v0, convergence_f converged,
            int max_iterations = 300, double deflation_tol = 1e-7,
            double delta = default_double, double eta = default_double,
            bool extended_local_reortho = true) {
  // Input:
  // mult:           multiplication function
  // v0:             starting vector for iterating
  // converged:      function to check for convergence
  // max_iterations: maximal number of iterations
  // deflation_tol:  tolerance for deflation
  // delta:          Desired level of orthogonality (default =
  //                 sqrt(eps/max_iterations)).
  // eta:            Level of orthogonality after reorthogonalization
  //                 (default = eps^(3/4)/sqrt(max_iterations)).
  // extended_local_reortho: flag whether to apply extended local
  // reorthogonalization
  using namespace arma;
  using vec_t = arma::Col<coeff_t>;
  using mat_t = arma::Mat<coeff_t>;

  lanczos_result<coeff_t> result;
  result.num_iterations = max_iterations;
  result.criterion = "unknown";

  // Set default orthogonality levels
  double eps = std::numeric_limits<double>::epsilon();
  if (delta == default_double) {
    delta = sqrt(eps / max_iterations);
  }
  if (eta == default_double) {
    eta = std::pow(eps, 0.75) / sqrt(max_iterations);
  }

  int64_t N = v0.size();
  double eps1 = sqrt((double)N) * eps / 2.0;
  double gamma = 1 / sqrt(2);

  vec &alphas = result.alphas;
  vec &betas = result.betas;
  alphas = zeros(max_iterations + 1);
  betas = zeros(max_iterations + 1);

  // Allocate Lanczos vector matric
  mat_t &V = result.V;
  try {
    V = zeros(N, max_iterations);
  } catch (...) {
    Log.err("Error in lanczos_pro: Unable to allocate Lanczos vector matrix, "
            "dim=({},{})",
            N, ax_iterations);
  }

  // Allocate Lanczos iteration vectors
  vec_t q, q_old;
  vec_t w;
  try {
    // Do I need all of them?
    q = zeros(N);
    q_old = zeros(N);
    w = zeros(N);
  } catch (...) {
    Log.err(
        "Error in lanczos_pro: Unable to allocate Lanczos iteration vectors, "
        "3 x dim={}",
        N);
  }

  // vectors to keep track of orthogonality
  vec &omegas = result.omegas;
  omegas = zeros(max_iterations);
  vec omegas_max = omegas;
  vec omegas_old = omegas;

  // Initialize with starting vector v0
  vec_t &r = v0;

  // Lanczos iterations
  for (int iteration = 0; iterations < max_iterations; ++iteration) {
    q_old = q;
    if (abs(betas(iteration)) < deflation_tol) {
      Log(1, "deflation detected at iteration {}", iteration);
      q = r;
      result.criterion = "deflation";
      break;
    } else {
      q = r / betas(iteration);
    }
    V.col(iteration) = q;

    // MVM
    mult(q, w);

    // orthogonalization
    r = w - beta(iteration) * q_old;
    alpha(iteration) = cdot(q, r);
    r = r - alpha(iteration) * q;

    // Extended local reorthogonalization (is this useful?)
    beta(iteration + 1) = sqrt(cdot(r, r));
    if (extended_local_reortho &&
        (beta(iteration + 1) < gamma * beta(iteration))) {
      if (iteration == 0) {
        double t1 = 0.;
        for (int i = 0; i < 2; ++i) {
          double t = cdot(q, r);
          r = r - q * t;
          t1 = t1 + t;
        }
        alpha(iteration) = alpha(iteration) + t1;
      } else {
        double t1 = cdot(q_old, r);
        double t2 = cdot(q, r);
        r = r - (q_old * t1 + q * t2);

        if (abs(beta(iteration)) > deflation_tol) {
          beta(iteration) += t1;
        }
        alpha(iteration) += t2;
      }
    }

    // Compute estimate for operator 2-norm (could be updated)
    mat Tmat = diagmat(alpha(span(0, iteration + 1))) +
               diagmat(betas(span(0, iteration)), -1) +
               diagmat(betas(span(0, iteration)), 1);
    double anorm = sqrt(norm(Tmat.trans() * Tmat, 1));

    // update omegas (orthogonalities), according to Simon algo
    // Horst D. Simon, The Lanczos algorithm with partial reorthogonalization
    if ((iteration > 0) && abs(beta(iteration + 1)) > deflation_tol) {
      double T = eps1 * anorm;
      double binv = 1 / beta(iteration + 1);
      omegas_old = omegas;
      omegas_old(0) = betas(1) * omegas(1) +
                      (alphas(0) - alphas(iteration)) * omegas(0) -
                      betas(iteration) * omegas_old(0);
      omegas_old(0) = binv * (omegas_old(0) + sign(omegas_old(0)) * T);
      for (int k = 1; k < iteration - 2; ++k) {
        omegas_old(k) = betas(k + 1) * omegas(k + 1) +
                        (alphas(k) - alphas(iteration)) * omegas(k) +
                        betas(k) * omegas(k - 1) -
                        betas(iteration, omega_old(k));
      }
      for (int k = 1; k < iteration - 2; ++k) {
        omegas_old(k) = binv * (omegas_old(k) + sign(omegas_old(k)) * T);
      }
      omegas_old(iteration - 1) = binv * T;
      swap(omegas, omegas_old);
      omegas(iteration) = eps1;
    }
  } // Lanczos iterations

  mat Tmat = diagmat(alpha(span(0, iteration + 1))) +
             diagmat(betas(span(0, iteration)), -1) +
             diagmat(betas(span(0, iteration)), 1);
  eig_sym(result.eigenvalues, Tmat);
  result.residual = r;
  
  return result;
}

} // namespace hydra
