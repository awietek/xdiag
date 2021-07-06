#pragma once

#include <lila/all.h>
#include <hydra/linalg/lanczos/tmatrix.h>

namespace hydra {

template <class coeff_t>
bool ConvergedEigenvalues(Tmatrix<coeff_t> const &tmat, int n_eigenvalue,
                          lila::real_t<coeff_t> precision);

} // namespace hydra
