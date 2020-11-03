#ifndef HYDRA_ENTANGLEMENT_ENTANGLEMENT_ENTROPY_H_
#define HYDRA_ENTANGLEMENT_ENTANGLEMENT_ENTROPY_H_

#include <lila/all.h>
#include <hydra/models/tjmodel.h>

namespace hydra {

template <class coeff_t>
lila::real_t<coeff_t> EntanglementEntropy(TJModel<coeff_t> const &model,
					  lila::Vector<coeff_t> const &wf,
					  int n_A_sites);

} // namespace hydra

#endif
