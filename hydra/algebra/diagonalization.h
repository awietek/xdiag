#pragma once

#include <lila/all.h>
#include <tuple>
#include <vector>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class Block>
lila::Vector<double>
LanczosEigenvalues_real(BondList const &bondlist, Couplings const &couplings,
                        Block &&block, int num_eigenvalue = 1,
                        double precision = 1e-12, int seed = 42,
                        int max_iterations = 1000) {
  using namespace lila;
  using vector_t = lila::Vector<double>;
  auto multiply = [&bondlist, &couplings, &block](vector_t const &v,
                                                  vector_t &w) {
    apply(bondlist, couplings, block, v, block, w);
  };

  // Initialize random starting vector
  idx_t dim = block.size();
  auto startstate = vector_t(dim);
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, seed);
  Random(startstate, gen, true);
  Normalize(startstate);

  // Run Lanczos
  auto res = LanczosEigenvalues(multiply, startstate, precision, num_eigenvalue,
                                "Ritz", max_iterations);
  return res.eigenvalues;
}

template <class Block>
lila::Vector<double>
LanczosEigenvalues_cplx(BondList const &bondlist, Couplings const &couplings,
                        Block &&block, int num_eigenvalue = 1,
                        double precision = 1e-12, int seed = 42,
                        int max_iterations = 1000) {
  using namespace lila;
  using vector_t = lila::Vector<complex>;
  auto multiply = [&bondlist, &couplings, &block](vector_t const &v,
                                                  vector_t &w) {
    apply(bondlist, couplings, block, v, block, w);
  };

  // Initialize random starting vector
  idx_t dim = block.size();
  auto startstate = vector_t(dim);
  normal_dist_t<complex> dist(0., 1.);
  normal_gen_t<complex> gen(dist, seed);
  Random(startstate, gen, true);
  Normalize(startstate);

  // Run Lanczos
  auto res = LanczosEigenvalues(multiply, startstate, precision, num_eigenvalue,
                                "Ritz", max_iterations);
  return res.eigenvalues;
}

template <class Block>
double e0_real(BondList const &bondlist, Couplings const &couplings,
               Block &&block, double precision = 1e-12, int seed = 42,
               int max_iterations = 1000) {
  auto eigs = LanczosEigenvalues_real(bondlist, couplings, block, 1, precision,
                                      seed, max_iterations);
  return eigs(0);
}

template <class Block>
double e0_cplx(BondList const &bondlist, Couplings const &couplings,
               Block &&block, double precision = 1e-12, int seed = 42,
               int max_iterations = 1000) {
  auto eigs = LanczosEigenvalues_cplx(bondlist, couplings, block, 1, precision,
                                      seed, max_iterations);
  return eigs(0);
}

} // namespace hydra
