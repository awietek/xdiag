#pragma once

#include <hydra/common.h>
#include <hydra/indexing/indexing_symmetric.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <class bit_t = std_bit_t,
          class GroupAction = PermutationGroupLookup<bit_t>>
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
  inline Representation const &irrep() const { return irrep_; }

  inline idx_t size() const { return indexing_.size(); }

  bool operator==(SpinhalfSymmetric const &rhs) const;
  bool operator!=(SpinhalfSymmetric const &rhs) const;

  // private:
  int n_sites_;
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;

  PermutationGroup permutation_group_;
  Representation irrep_;

  indexing::IndexingSymmetric<bit_t, GroupAction> indexing_;
};

} // namespace hydra
