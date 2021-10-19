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

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indexspinhalf.h>
#include <hydra/indexing/indexsymmetrized.h>
#include <hydra/utils/bitops.h>
#include <hydra/utils/range.h>

#include "spinlessfermions.h"

namespace hydra {

SpinlessFermions::SpinlessFermions(
    int n_sites, const std::vector<std::pair<int, int>> hoppings,
    const std::vector<std::pair<int, int>> interactions)
    : n_sites_(n_sites), hoppings_(hoppings), interactions_(interactions) {}

std::vector<int> SpinlessFermions::quantumnumbers() {
  std::vector<int> qns;
  for (int np = 0; np <= n_sites_; ++np)
    qns.push_back(np);
  return qns;
}

lila::Matrix<double> SpinlessFermions::matrix(double t, double V,
                                              int qn) const {
  using hilbertspaces::Spinhalf;
  using indexing::IndexSpinhalf;
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  Spinhalf<uint64> hs(n_sites_, qn);
  IndexSpinhalf<uint64, uint64> indexing(hs);
  int dim = indexing.size();
  lila::Matrix<double> hamilton(dim, dim);
  lila::Zeros(hamilton);

  // Apply hopping
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
    int idx = 0;
    for (auto state : hs) {
      if (((state & flipmask) != 0) && ((state & flipmask) != flipmask)) {
        const double fermi =
            popcnt(gbits(state, s2 - s1 - 1, s1 + 1)) % 2 == 0 ? 1. : -1.;
        auto new_state = state ^ flipmask;
        int new_idx = indexing.index(new_state);
        hamilton(new_idx, idx) -= fermi * t;
      }
      ++idx;
    }
  }

  // Apply interaction
  for (auto pair : interactions_) {
    int s1 = pair.first;
    int s2 = pair.second;
    int idx = 0;
    for (auto state : hs) {
      hamilton(idx, idx) +=
          V * ((double)gbit(state, s1)) * ((double)gbit(state, s2));
      ++idx;
    }
  }

  return hamilton;
}

lila::Matrix<complex>
SpinlessFermions::matrix(double t, double V, int qn,
                         CharacterTable &character_table,
                         std::string representation_name) const {
  using hilbertspaces::Spinhalf;
  using indexing::IndexSymmetrized;
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  Spinhalf<uint64> hs(n_sites_, qn);
  IndexSymmetrized<Spinhalf<uint64>, uint32, true> indexing(
      hs, character_table, representation_name);
  std::vector<complex> characters =
      character_table.characters(representation_name);
  int dim = indexing.size();
  lila::Matrix<complex> hamilton(dim, dim);
  lila::Zeros(hamilton);

  auto space_group = character_table.little_group(representation_name);

  // Construct hopping matrix
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);

    uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);

      // Apply hopping terms
      if (((state & flipmask) != 0) && ((state & flipmask) != flipmask)) {
        auto new_state = state ^ flipmask;
        auto representative = new_state;
        int n_sym = indexing.find_representative(&representative);
        int new_idx = indexing.index(representative);

        if (new_idx != -1) {
          double fermi =
              popcnt(gbits(state, s2 - s1 - 1, s1 + 1)) % 2 == 0 ? 1. : -1.;
          double fermi2 =
              space_group.fermi_sign<Spinhalf<uint64>>(n_sym, new_state);
          complex coeff = -t * fermi * fermi2 * characters[n_sym] *
                          indexing.norm(new_idx) / indexing.norm(idx);
          hamilton(new_idx, idx) += coeff;
        }
      }

    } // loop over basis states
  }   // loop over hoppings

  // Construct interaction matrix
  for (auto pair : interactions_) {
    int s1 = pair.first;
    int s2 = pair.second;
    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);
      hamilton(idx, idx) +=
          V * ((double)gbit(state, s1)) * ((double)gbit(state, s2));
    }
  }
  return hamilton;
}

} // namespace hydra
