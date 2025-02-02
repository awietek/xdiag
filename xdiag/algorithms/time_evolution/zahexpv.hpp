#pragma once
//
// Created by Luke Staszewski on 27.01.23.
//

#include <cmath>
#include <tuple>
#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/algorithms/norm_estimate.hpp>
#include <xdiag/algorithms/time_evolution/expm.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>

#include <xdiag/parallel/mpi/allreduce.hpp>

namespace xdiag {

template <typename apply_A_f, typename dot_f>
inline std::tuple<double, double>
zahexpv(double time, apply_A_f &&apply_A, dot_f &&dot, arma::cx_vec &w,
        double anorm, double tol = 1e-12, int m = 30) try {

  /* perform the time evolution via mat exponential applied to vector as
  outlined in expokit paper i.e. returns w = eˆ(At)*v, where A is anti-hermitian
  matrix and t is a real parameter params: time: real, amount of time to time
  evolve initial state v by apply_A: function for performing mat multiplication
  on-the-fly, apply_A(v): v -> w = Av v: initial vector to be time evolved
      anorm: an estimate for the norm of H outlined in armadillo docs for
  norm(X, p) p = 'inf'
  krylov space - can be updated but will set to 30 or less probably
  */

  auto norm = [&dot](arma::cx_vec const &v) {
    return std::sqrt(xdiag::real(dot(v, v)));
  };

  //     parameter initialisation
  int64_t n = w.size();
  //    m = std::min(m, n);
  int mxrej = 10;
  double btol = 1.0e-7;
  double gamma = 0.9;
  double delta = 1.2; // recommended but can be adjusted
  int mb = m;
  double t_out = std::abs(time);
  double nstep = 0;
  double t_new = 0;
  double t_now = 0;
  double s_error = 0;

  // getting machine precision could just use built in
  double p1 = 4.0 / 3.0;
  double eps = 0;
  while (eps == 0) {
    double p2 = p1 - 1.0;
    double p3 = p2 + p2 + p2;
    eps = std::abs(p3 - 1.0);
  }
  if (tol < eps)
    tol = sqrt(eps);
  double rndoff = eps * anorm;

  // determining the first time step
  int k1 = 2;
  double xm = 1 / (double)m;
  double normv = norm(w);
  double beta = normv;
  double fact = std::pow((m + 1) / exp(1),
                         (m + 1)) *
                sqrt(2 * M_PI * (m + 1)); // sterling's approx.
  t_new = (1 / anorm) * std::pow((fact * tol) / (4 * beta * anorm), xm);

  double s = std::pow(10, floor(log10(t_new)) - 1);
  t_new = ceil(t_new / s) * s;
  int sgn = arma::sign(time);
  nstep = 0;

  // hump determines if the matrix is conditioned or not ( <1 => well
  // conditioned ) expokit paper still stipulates algorithm can work even if if
  // hump > 1
  double hump = normv;

  // the actual time evolution
  while (t_now < t_out) {
    nstep++;

    double t_step = std::min(t_out - t_now, t_new);

    // storage of the basis vectors for the Krylov space
    arma::cx_mat V(n, m + 1, arma::fill::zeros);

    // A written in the basis of the krylov space
    arma::mat H(m + 2, m + 2, arma::fill::zeros);

    V(arma::span(0, n - 1), 0) = (1 / beta) * w;

    // building up the basis vectors of the Krylov space (the lanczos steps)
    //  n.b. this is done only once for each time step
    arma::cx_vec p;
    for (int j = 0; j < m; j++) {
      p = apply_A(
          arma::conv_to<arma::cx_vec>::from(V(arma::span(0, n - 1), j)));

      H(j, j) = real(dot(V(arma::span(0, n - 1), j), p));
      p -= H(j, j) * V(arma::span(0, n - 1), j);
      if (j != 0) {
        H(j - 1, j) =
            -H(j, j - 1); // this step only works for anti- hermitian mat
        p -= H(j - 1, j) * V(arma::span(0, n - 1), j - 1);
      }
      s = norm(p);

      if (s <= btol) { // happy breakdown
        Log(2, "time evolution using zahexpv: happy breakdown for j = {}", j);
        k1 = 0;
        mb = j + 1;
        t_step = t_out - t_now;
        break;
      }
      H(j + 1, j) = s;
      V(arma::span(0, n - 1), j + 1) = (1 / s) * p;
    } // end of loop for building up krylov space basis vectors

    double avnorm = 0;
    if (k1 != 0) { // i.e. no happy break-down
      H(m + 1, m) = 1;
      avnorm = norm(apply_A(
          arma::conv_to<arma::cx_vec>::from(V(arma::span(0, n - 1), m))));
    }
    // iteratively selecting the step size until desired tol achieved
    int ireject = 0; // won't do more than 10 iterations per time step
    int mx = mb + k1;
    arma::mat F;
    double err_loc = 1;

    while (ireject <= mxrej) {

      F = expm(arma::mat(sgn * t_step * H.submat(0, 0, mx - 1, mx - 1)));

      if (k1 == 0) {    // case of happy break-down
        err_loc = btol; // TODO: try tol vs btol :this was error_local in matlab
                        // imp. should it be err_lco
        break;
      }

      double phi1 = std::abs(beta * F(m, 0));
      double phi2 = std::abs(beta * F(m + 1, 0) * avnorm);

      if (phi1 > 10 * phi2) {
        err_loc = phi2;
        xm = 1 / (double)m;
      } else if (phi1 > phi2) {
        err_loc = (phi1 * phi2) / (phi1 - phi2);
        xm = 1 / (double)m;
      } else {
        err_loc = phi1;
        xm = 1 / (double)(m - 1);
      }

      if (err_loc <= delta * t_step * tol) {
        break;
      }

      t_step = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
      s = std::pow(10, (floor(log10(t_step)) - 1));
      t_step = ceil(t_step / s) * s;
      if (ireject == mxrej) {
        XDIAG_THROW(
            "The requested tolerance is too high (irej > "
            "10). Try increasing m (Krylov space dimension) or decreasing "
            "the precision.");
      }
      ireject++;

    } // end of while loop for checking whether time step has converged
      // if happy - breakdown size of H and F taken as m, and m+1 otherwise

    mx = mb + std::max(0, k1 - 1);
    // get first column of F := F0
    arma::vec F0 = F(arma::span(0, mx - 1), 0);
    w.zeros();
    w = beta * V(arma::span(0, n - 1), arma::span(0, mx - 1)) * F0;
    beta = norm(w);
    hump = std::max(hump, beta);

    t_now += t_step;

    t_new = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
    s = std::pow(10, floor(log10(t_new)) - 1);
    t_new = ceil(t_new / s) * s;

    err_loc = std::max(err_loc, rndoff);
    s_error = s_error + err_loc;
  } // end of full time evolution

  double err = s_error;
  hump = hump / normv;
  Log(1, "zaexph finished: # steps = {}, # MVM = {}, est. error: {}, hump: {}",
      nstep, nstep * m, err, hump);
  return {err, hump};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return {0., 0.};
}

} // namespace xdiag
