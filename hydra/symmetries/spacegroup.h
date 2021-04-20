#pragma once

#include <string>
#include <vector>

#include <hydra/symmetries/symmetry_operations.h>

namespace hydra {

class SpaceGroup {
public:
  SpaceGroup(std::vector<std::vector<int>> const &symmetries);

  template <class bit_t> inline bit_t apply(int n_sym, bit_t state) const {
    return detail::apply_permutation(
        state, n_sites_, symmetries_internal_.data() + n_sym * n_sites_);
  }

  template <class bit_t>
  inline double fermi_sign(int n_sym, bit_t state) const {
    return detail::fermi_sign(state, n_sites_,
                              symmetries_internal_.data() + n_sym * n_sites_);
  }

  SpaceGroup subgroup(std::vector<int> const &symmetry_numbers) const;

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

SpaceGroup read_spacegroup(std::string filename);

} // namespace hydra
