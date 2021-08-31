#pragma once

#include <unordered_map>

#include <hydra/common.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <class bit_t = std_bit_t, class GroupAction = PermutationGroupAction>
class SpinhalfSymmetric {
public:
  SpinhalfSymmetric() = default;
  SpinhalfSymmetric(int n_sites, int n_up, PermutationGroup permutation_group,
                    Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_up_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline GroupAction const &group_action() const { return group_action_; }
  inline Representation const &irrep() const { return irrep_; }

  inline idx_t size() const { return size_; }

  bool operator==(SpinhalfSymmetric const &rhs) const;
  bool operator!=(SpinhalfSymmetric const &rhs) const;

// private:
  idx_t index(bit_t state) const;

  int n_sites_;

  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;

  PermutationGroup permutation_group_;
  GroupAction group_action_;
  Representation irrep_;

  std::vector<bit_t> states_;
  std::vector<double> norms_;
  idx_t size_;
};

} // namespace hydra
