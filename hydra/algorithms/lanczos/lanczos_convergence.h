#pragma once

#include "extern/armadillo/armadillo"
#include <hydra/algorithms/lanczos/tmatrix.h>

namespace hydra::lanczos {

bool converged_eigenvalues(Tmatrix const &tmat, int n_eigenvalue,
                           double precision);
bool converged_time_evolution(Tmatrix const &tmat, complex tau,
                              double precision, double nrm);

} // namespace hydra
