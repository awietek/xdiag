#include "entanglement_entropy.h"

#include <cmath>
#include <hydra/entanglement/reduced_density_matrix.h>

namespace hydra {

template <class coeff_t>
lila::real_t<coeff_t> EntanglementEntropy(TJModel<coeff_t> const &model,
                                          lila::Vector<coeff_t> const &wf,
                                          int n_A_sites) {
  using model_t = TJModel<coeff_t>;
  using qn_t = typename model_t::qn_t;

  assert(n_A_sites <= model.n_sites());

  auto qn_tot = model.qn();
  auto n_B_sites = model.n_sites() - n_A_sites;

  lila::real_t<coeff_t> entropy = 0;
  for (number_t n_A_up = 0; n_A_up <= n_A_sites; ++n_A_up)
    for (number_t n_A_dn = 0; n_A_dn <= n_A_sites; ++n_A_dn) {
      auto qn_A = qn_t({n_A_up, n_A_dn});
      auto qn_B = qn_tot - qn_A;

      if (valid(qn_A, n_A_sites) && valid(qn_B, n_B_sites)) {
        // std::cout << "Tot  " << String(qn_tot) << std::endl;
        // std::cout << "A    " << String(qn_A) << std::endl;
        // std::cout << "B    " << String(qn_B) << std::endl;

        auto rho_A = ReducedDensityMatrix(model, wf, qn_A, n_A_sites);
        assert(lila::close(rho_A, lila::Herm(rho_A)));

        auto rho_eigs = lila::EigenvaluesSym(rho_A);
        for (auto eig : rho_eigs) {
          if (!lila::close(eig, 0.))
            entropy -= eig * log(eig);
        }
      }
    }
  return entropy;
}

template double EntanglementEntropy(TJModel<double> const &model,
                                    lila::Vector<double> const &wf,
                                    int n_A_sites);

template double EntanglementEntropy(TJModel<complex> const &model,
                                    lila::Vector<complex> const &wf,
                                    int n_A_sites);
} // namespace hydra
