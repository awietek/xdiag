#pragma once

#include <xdiag/common.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

struct EigsLanczosResult {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  State eigenvectors;
  int64_t niterations;
  std::string criterion;
};

XDIAG_API EigsLanczosResult eigs_lanczos(OpSum const &ops, Block const &block,
                                         int64_t neigvals = 1,
                                         double precision = 1e-12,
                                         int64_t max_iterations = 1000,
                                         double deflation_tol = 1e-7,
                                         int64_t random_seed = 42);

XDIAG_API EigsLanczosResult eigs_lanczos(OpSum const &ops, State const &state0,
                                         int64_t neigvals = 1,
                                         double precision = 1e-12,
                                         int64_t max_iterations = 1000,
                                         double deflation_tol = 1e-7);

} // namespace xdiag
