#ifndef HYDRA_ENTANGLEMENT_REDUCED_DENSITY_MATRIX_H_
#define HYDRA_ENTANGLEMENT_REDUCED_DENSITY_MATRIX_H_

#include <lila/all.h>
#include <hydra/models/tjmodel.h>

namespace hydra {

template <class coeff_t>
lila::Matrix<coeff_t> ReducedDensityMatrix(TJModel<coeff_t> const &model,
                                           lila::Vector<coeff_t> const &wf,
					   qn_tj qn, int n_A_sites);

} // namespace hydra

#endif
