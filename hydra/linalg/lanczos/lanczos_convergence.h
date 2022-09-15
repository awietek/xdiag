#pragma once

#include "extern/armadillo/armadillo"
#include <hydra/linalg/lanczos/tmatrix.h>

namespace hydra {

bool ConvergedEigenvalues(Tmatrix const &tmat, int n_eigenvalue,
                          double precision);

} // namespace hydra
