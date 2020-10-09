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

#ifndef HYDRA_SYMMETRIES_SPACEGROUP_
#define HYDRA_SYMMETRIES_SPACEGROUP_

#include <string>
#include <vector>

#include <hydra/symmetries/symmetrydetail.h>

namespace hydra {

class SpaceGroup {
public:
  SpaceGroup(const std::vector<std::vector<int>> &symmetries);

  template <class state_t>
  inline state_t apply(int n_sym, state_t const &state) const {
    return detail::apply_permutation(
        state, n_sites_, symmetries_internal_.data() + n_sym * n_sites_);
  }

  template <class hilbertspace_t>
  inline double
  fermi_sign(const int &n_sym,
             const typename hilbertspace_t::state_t &state) const {
    return detail::fermi_sign<hilbertspace_t>(
        state, n_sites_, symmetries_internal_.data() + n_sym * n_sites_);
  }

  SpaceGroup subgroup(const std::vector<int> &symmetry_numbers) const;

  int n_sites() const { return n_sites_; }
  int n_symmetries() const { return n_symmetries_; }
  const std::vector<std::vector<int>> &symmetries() const {
    return symmetries_;
  }

private:
  int n_sites_;
  int n_symmetries_;
  std::vector<std::vector<int>> symmetries_;
  std::vector<int> symmetries_internal_; // size = n_symmetries_*n_sites_
};

// void Print(const SpaceGroup& group);
SpaceGroup read_spacegroup(std::string filename);

} // namespace hydra

#endif
