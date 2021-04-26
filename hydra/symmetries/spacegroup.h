#pragma once

#include <string>
#include <vector>

#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/spacegroup_operator.h>

namespace hydra {

template <class bit_t = std_bit_t,
          class SpaceGroupOperator = SpaceGroupOperator<bit_t>>
class SpaceGroup {
public:
  SpaceGroup() = default;
  explicit SpaceGroup(std::vector<std::vector<int>> const &symmetries);

  inline bit_t apply(int sym, bit_t state) const {
    return spacegroup_operator_.apply(sym, state);
  }
  inline bit_t representative(bit_t state) const {
    return spacegroup_operator_.representative(state);
  }
  inline std::tuple<bit_t, int> representative_index(bit_t state) const {
    return spacegroup_operator_.representative_index(state);
  }
  inline std::tuple<bit_t, int, const int *>
  representative_indices(bit_t state) const {
    return spacegroup_operator_.representative_indices(state);
  }
  inline double fermi_sign(int sym, bit_t state) const {
    return spacegroup_operator_.fermi_sign(sym, state);
  }

  SpaceGroup<bit_t, SpaceGroupOperator>
  subgroup(std::vector<int> const &symmetry_numbers) const;

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline int permutation(int sym, int site) const {
    return permutation_array_[n_sites_ * sym + site];
  }

  bool operator==(SpaceGroup const& rhs) const;
  bool operator!=(SpaceGroup const& rhs) const;
  
private:
  int n_sites_;
  int n_symmetries_;
  SpaceGroupOperator spacegroup_operator_;
  std::vector<std::vector<int>> symmetries_;
  std::vector<int> permutation_array_; // size = n_symmetries_*n_sites_
};

} // namespace hydra
