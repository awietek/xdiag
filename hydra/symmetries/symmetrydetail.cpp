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

#include <algorithm>
#include <vector>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/utils/bitops.h>

#include "symmetrydetail.h"

namespace hydra {
namespace detail {

bool is_valid_permutation(const std::vector<int> &permutation) {
  for (int i = 0; i < (int)permutation.size(); ++i)
    if (std::find(permutation.begin(), permutation.end(), i) ==
        permutation.end())
      return false;
  return true;
}

template <class state_t>
double fermi_sign_spinhalf(const state_t &state, const int &n_sites,
                           const int *permutation) {
  std::vector<int> sort(n_sites);
  int sum = 0;
  int nf, k;

  for (nf = 0, k = 0; k < n_sites; ++k)
    if (utils::gbit(state, k)) {
      sort[nf] = permutation[k];
      ++nf;
    }

  int old_sum;
  while (true) {
    old_sum = sum;
    for (k = 0; k < (nf - 1); ++k)
      if (sort[k + 1] < sort[k]) {
        sum++;
        std::swap(sort[k + 1], sort[k]);
      }
    if (old_sum == sum)
      break;
  }
  return (sum % 2 ? -1 : 1);
}

using hilbertspaces::Spinhalf;

template <>
double fermi_sign<Spinhalf<uint16>>(const uint16 &state, const int &n_sites,
                                    const int *permutation) {
  return fermi_sign_spinhalf(state, n_sites, permutation);
}

template <>
double fermi_sign<Spinhalf<uint32>>(const uint32 &state, const int &n_sites,
                                    const int *permutation) {
  return fermi_sign_spinhalf(state, n_sites, permutation);
}

template <>
double fermi_sign<Spinhalf<uint64>>(const uint64 &state, const int &n_sites,
                                    const int *permutation) {
  return fermi_sign_spinhalf(state, n_sites, permutation);
}
} // namespace detail
} // namespace hydra
