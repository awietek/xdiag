#ifndef HYDRA_ENTANGLEMENT_ENTANGLEMENT_ENTROPY_H_
#define HYDRA_ENTANGLEMENT_ENTANGLEMENT_ENTROPY_H_

#include <hydra/models/tjmodel.h>
#include <lila/all.h>

namespace hydra {

template <class model_t>
lila::real_t<typename model_t::coeff_t>
EntanglementEntropy(model_t const &model, typename model_t::vector_t const &wf,
                    int n_A_sites);

} // namespace hydra

#endif
