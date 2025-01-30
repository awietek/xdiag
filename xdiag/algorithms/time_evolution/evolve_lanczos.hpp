#pragma once

#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag {

struct evolve_lanczos_result_t {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
  State state;
};

XDIAG_API evolve_lanczos_result_t
evolve_lanczos(OpSum const &H, State psi, double tau, double precision = 1e-12,
               double shift = 0., bool normalize = false,
               int64_t max_iterations = 1000, double deflation_tol = 1e-7);

XDIAG_API evolve_lanczos_result_t
evolve_lanczos(OpSum const &H, State psi, complex tau, double precision = 1e-12,
               double shift = 0., bool normalize = false,
               int64_t max_iterations = 1000, double deflation_tol = 1e-7);

struct evolve_lanczos_inplace_result_t {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
};

XDIAG_API evolve_lanczos_inplace_result_t evolve_lanczos_inplace(
    OpSum const &H, State &psi, double tau, double precision = 1e-12,
    double shift = 0., bool normalize = false, int64_t max_iterations = 1000,
    double deflation_tol = 1e-7);

XDIAG_API evolve_lanczos_inplace_result_t evolve_lanczos_inplace(
    OpSum const &H, State &psi, complex tau, double precision = 1e-12,
    double shift = 0., bool normalize = false, int64_t max_iterations = 1000,
    double deflation_tol = 1e-7);

} // namespace xdiag
