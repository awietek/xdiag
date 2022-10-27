#pragma once

#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/blocks/blocks.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

// Lanczos which  overwrites starting vector v0
template <class coeff_t>
Tmatrix
lanczos_eigenvalues_inplace(BondList const &bonds, State<coeff_t> &state_0,
                            int num_eigenvalue = 0, double precision = 1e-12,
                            int max_iterations = 1000,
                            double deflation_tol = 1e-7);

// Lanczos which does not overwrite v0
template <class coeff_t>
Tmatrix lanczos_eigenvalues(BondList const &bonds, State<coeff_t> state_0,
                            int num_eigenvalue = 0, double precision = 1e-12,
                            int max_iterations = 1000,
                            double deflation_tol = 1e-7);

// Lanczos which does not overwrite v0
template <class coeff_t>
Tmatrix lanczos_eigenvalues(BondList const &bonds, Block const &block,
                            int num_eigenvalue = 0, double precision = 1e-12,
                            int max_iterations = 1000,
                            double deflation_tol = 1e-7);

Tmatrix lanczos_eigenvalues_real(BondList const &bonds, Block const &block,
                                 int num_eigenvalue = 0,
                                 double precision = 1e-12, uint64_t seed = 42,
                                 int max_iterations = 1000,
                                 double deflation_tol = 1e-7);

Tmatrix lanczos_eigenvalues_cplx(BondList const &bonds, Block const &block,
                                 int num_eigenvalue = 0,
                                 double precision = 1e-12, uint64_t seed = 42,
                                 int max_iterations = 1000,
                                 double deflation_tol = 1e-7);

} // namespace hydra
