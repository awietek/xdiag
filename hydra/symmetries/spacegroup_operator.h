#pragma once

#include <tuple>
#include <vector>

#include <hydra/common.h>
#include <hydra/symmetries/spacegroup_operations.h>

namespace hydra {

template <class bit_t> class SpaceGroupOperator {
public:
  SpaceGroupOperator() = default;
  SpaceGroupOperator(int n_sites, std::vector<int> const &permutation_array);

  inline bit_t apply(int sym, bit_t state) const {
    return detail::apply_permutation(
        state, n_sites_, permutation_array_.data() + sym*n_sites_);
  }

  bit_t representative(bit_t state) const;
  std::tuple<bit_t, int> representative_index(bit_t state) const;

private:
  int n_sites_;
  int n_sym_;
  std::vector<int> permutation_array_;
};

} // namespace hydra
