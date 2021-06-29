#ifndef HYDRA_ENTANGLEMENT_REDUCED_DENSITY_MATRIX_H_
#define HYDRA_ENTANGLEMENT_REDUCED_DENSITY_MATRIX_H_

#include <hydra/models/tjmodel.h>
#include <hydra/models/tjmodelmpi.h>
#include <lila/all.h>

namespace hydra {

template <class coeff_t>
lila::Matrix<coeff_t> ReducedDensityMatrix(TJModel<coeff_t> const &model,
                                           int n_A_sites, qn_tj qn,
                                           lila::Vector<coeff_t> const &wf);

template <class coeff_t>
lila::Matrix<coeff_t> ReducedDensityMatrix(TJModelMPI<coeff_t> const &model,
					   int n_A_sites, qn_tj qn,
					   lila::VectorMPI<coeff_t> const &wf);

} // namespace hydra

#endif
