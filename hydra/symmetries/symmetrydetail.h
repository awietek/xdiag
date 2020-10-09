// Copyright 2018 Alexander Wietek - All Rights Reserved.
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

#ifndef HYDRA_SYMMETRIES_SYMMETRYDETAIL_
#define HYDRA_SYMMETRIES_SYMMETRYDETAIL_

#include <vector>

namespace hydra {
namespace detail {

bool is_valid_permutation(const std::vector<int> &permutation);

template <class state_t>
state_t apply_permutation(state_t const &state, int const &n_sites,
                          const int *permutation) {
  state_t tstate = 0;
  for (int i = 0; i < n_sites; ++i) {
    const int val = siteval(state, i);
    set_siteval(&tstate, permutation[i], val);
  }
  return tstate;
}

template <class hilbertspace_t>
double fermi_sign(const typename hilbertspace_t::state_t &state,
                  const int &n_sites, const int *permutation);

} // namespace detail
} // namespace hydra

#endif
