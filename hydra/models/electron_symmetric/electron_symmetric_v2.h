#pragma once

#include <map>

#include <hydra/common.h>
#include <hydra/indexing/indexing_symmetric_fermionic.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <class bit_t = std_bit_t, class GroupAction = PermutationGroupLookup>
class ElectronSymmetricV2 {
public:
  ElectronSymmetricV2() = default;
  ElectronSymmetricV2(int n_sites, int nup, int ndn,
                      PermutationGroup permutation_group,
		      Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_up_; }
  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline GroupAction const &group_action() const { return group_action_; }
  inline Representation const &irrep() const { return irrep_; }

  inline idx_t size() const { return size_; }

  bool operator==(ElectronSymmetricV2 const &rhs) const;
  bool operator!=(ElectronSymmetricV2 const &rhs) const;

private:
  int n_sites_;

  bool charge_conserved_;
  int charge_;
  bool sz_conserved_;
  int sz_;
  int n_up_;
  int n_dn_;

  PermutationGroup permutation_group_;
  GroupAction group_action_;
  Representation irrep_;

  IndexingSymmetricFermionic<bit_t, GroupAction> idxes_up_;
  IndexingSymmetricFermionic<bit_t, GroupAction> idxes_dn_;

  std::vector<int> stabilizers_up_;
  std::vector<int> stabilizers_dn_;

  std::vector<std::pair<idx_t, idx_t>> stabilizers_up_limits_;
  std::vector<std::pair<idx_t, idx_t>> stabilizers_dn_limits_;

  std::vector<bit_t> dns_for_up_; 
  std::vector<bit_t> ups_for_dn_;  
  
  idx_t size_;
};

} // namespace hydra
