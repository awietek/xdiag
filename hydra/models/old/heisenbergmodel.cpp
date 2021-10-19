// Copyright 2019 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <lila/special.h>

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/indexing/index_spinhalf.h>
#include <hydra/indexing/index_symmetrized.h>

#include <hydra/symmetries/charactertable.h>

#include <hydra/utils/bitops.h>
#include <hydra/utils/range.h>

#include "heisenbergmodel.h"

namespace hydra {

template <class bit_t>
lila::Matrix<double> HeisenbergModel<bit_t>::matrix(BondList bondlist,
                                                    Couplings couplings,
                                                    int qn) const {
  using utils::gbit;
  using utils::popcnt;

  int n_sites = bondlist.n_sites();

  BasisSpinHalf<bit_t> hs(n_sites, qn);
  IndexSpinHalf<bit_t, bit_t> indexing(hs);
  int dim = indexing.size();
  lila::Matrix<double> hamilton(dim, dim);
  lila::Zeros(hamilton);

  auto heisenberg_bonds = bondlist.bonds_of_type("HEISENBERG");
  for (auto bond : heisenberg_bonds) {
    assert(bond.size() == 2); // Heisenberg bonds must have length 2
    int s1 = bond.sites(0);
    int s2 = bond.sites(1);
    std::string coupling = bond.coupling();
    double J = lila::real(couplings[coupling]);

    // Apply Heisenberg operator on sites s1, s2
    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    for (int idx = 0; idx < indexing.size(); ++idx) {
      auto state = indexing.state(idx);

      if (siteval(state, s1) == siteval(state, s2))
        hamilton(idx, idx) += J / 4.; // Ising
      else {
        hamilton(idx, idx) -= J / 4.; // Ising

        // Exchange term
        state_t new_state = {state.spins ^ flipmask};
        int new_idx = indexing.index(new_state);
        hamilton(new_idx, idx) += J / 2.;
      }
    }
  }
  return hamilton;
}

// template <class bit_t>
// lila::Matrix<complex>
// HeisenbergModel<bit_t>::matrix(BondList bondlist, Couplings couplings, int qn,
//                                CharacterTable &character_table,
//                                std::string representation_name) const {
//   using utils::gbit;
//   using utils::popcnt;
//   using utils::range;

//   int n_sites = bondlist.n_sites();

//   BasisSpinHalf<bit_t> hs(n_sites, qn);
//   IndexSymmetrized<BasisSpinHalf<bit_t>> indexing(hs, character_table,
//                                                   representation_name);
//   std::vector<complex> characters =
//       character_table.characters(representation_name);
//   int dim = indexing.size();
//   lila::Matrix<complex> hamilton(dim, dim);
//   lila::Zeros(hamilton);

//   auto heisenberg_bonds = bondlist.bonds_of_type("HEISENBERG");
//   for (auto bond : heisenberg_bonds) {
//     assert(bond.size() == 2); // Heisenberg bonds must have length 2
//     int s1 = bond.sites(0);
//     int s2 = bond.sites(1);
//     std::string coupling = bond.coupling();
//     double J = lila::real(couplings[coupling]);

//     // Apply Heisenberg operator on sites s1, s2
//     bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
//     for (int idx : range<>(indexing.size())) {
//       auto state = indexing.state(idx);

//       if (siteval(state, s1) == siteval(state, s2))
//         hamilton(idx, idx) += J / 4.; // Ising
//       else {
//         hamilton(idx, idx) -= J / 4.; // Ising

//         // Exchange term
//         state_t new_state = {state.spins ^ flipmask};
//         state_t representative = new_state;
//         int n_sym = indexing.find_representative(&representative);
//         int new_idx = indexing.index(representative);

//         if (new_idx != -1) {
//           complex character = characters[n_sym];
//           complex coeff =
//               indexing.norm(new_idx) / indexing.norm(idx) * character * J / 2.;
//           hamilton(new_idx, idx) += coeff;
//         }
//       }
//     }
//   }

//   return hamilton;
// }

// template class HeisenbergModel<uint16>;
template class HeisenbergModel<uint32>;
template class HeisenbergModel<uint64>;

} // namespace hydra
