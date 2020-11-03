#include "reduced_density_matrix.h"

#include <hydra/indexing/index_table.h>

namespace hydra {

template <class coeff_t>
lila::Matrix<coeff_t> ReducedDensityMatrix(TJModel<coeff_t> const &model,
                                           lila::Vector<coeff_t> const &wf,
                                           qn_tj qn_A, int n_sites_A) {
  using model_t = TJModel<coeff_t>;
  using idx_t = typename model_t::idx_t;
  using basis_t = typename model_t::basis_t;

  auto qn_tot = model.qn();
  auto qn_B = qn_tot - qn_A;
  int n_sites_tot = model.n_sites();
  int n_sites_B = n_sites_tot - n_sites_A;

  assert(n_sites_tot >= 0);
  assert(n_sites_A >= 0);
  assert(n_sites_B >= 0);

  assert(valid(qn_tot, n_sites_tot));
  assert(valid(qn_A, n_sites_A));
  assert(valid(qn_B, n_sites_B));

  auto basis_tot = basis_t(n_sites_tot, qn_tot);
  auto basis_A = basis_t(n_sites_A, qn_A);
  auto basis_B = basis_t(n_sites_B, qn_B);

  IndexTable<basis_t, idx_t> indexing(basis_tot);

  // Compute the partial trace over B
  auto dim = basis_A.size();
  auto rho_A = lila::Zeros<coeff_t>(dim, dim);

  idx_t i_A1 = 0;
  for (auto state_A1 : basis_A) {
    idx_t i_A2 = 0;
    for (auto state_A2 : basis_A) {
      for (auto state_B : basis_B) {
        auto state_B_shifted = (state_B << n_sites_A);
        auto state_1 = state_A1 | state_B_shifted;
        auto state_2 = state_A2 | state_B_shifted;
        auto idx_1 = indexing.index(state_1);
        auto idx_2 = indexing.index(state_2);

        // std::cout << i_A1 << " " << i_A2 << " " << idx_1 << " " << idx_2
        //           << "\n";
	// std::cout << String(state_A1, n_sites_A) << "\n";
	// std::cout << String(state_A2, n_sites_A) << "\n";
	// std::cout << String(state_B, n_sites_B) << "\n";
	
	// std::cout << String(state_1, n_sites_tot) << "\n";
	// std::cout << String(state_2, n_sites_tot) << "\n\n";
	
        rho_A(i_A1, i_A2) += lila::conj(wf(idx_1)) * wf(idx_2);
      }
      ++i_A2;
    }
    ++i_A1;
  }

  return rho_A;
}

template lila::Matrix<double>
ReducedDensityMatrix(TJModel<double> const &model,
                     lila::Vector<double> const &wf, qn_tj qn_A, int n_sites_A);

template lila::Matrix<complex>
ReducedDensityMatrix(TJModel<complex> const &model,
                     lila::Vector<complex> const &wf, qn_tj qn_A,
                     int n_sites_A);

} // namespace hydra
