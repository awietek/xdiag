#pragma once

#include <map>

#include <hydra/common.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <class bit_t = std_bit_t, class GroupAction = PermutationGroupAction>
class ElectronSymmetric {
public:
  ElectronSymmetric() = default;
  ElectronSymmetric(int n_sites, int nup, int ndn,
                    PermutationGroup permutation_group, Representation irrep);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_up_; }
  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline PermutationGroup const &permutation_group() const { return permutation_group_; }
  inline GroupAction const &group_action() const { return group_action_; }
  inline Representation const &irrep() const { return irrep_; }

  inline idx_t size() const { return size_; }
  inline bit_t up(idx_t idx) const { return ups_[idx]; }
  inline bit_t dn(idx_t idx) const { return dns_[idx]; }
  inline double norm(idx_t idx) const { return norms_[idx]; }

  std::tuple<bit_t, bit_t> representative(bit_t ups, bit_t dns) const;
  std::tuple<bit_t, bit_t, int> representative_index(bit_t ups,
                                                     bit_t dns) const;
  idx_t index(bit_t ups, bit_t dns) const;
  idx_t index_switch(bit_t ups, bit_t dns) const;
  inline idx_t index_switch_to_index(idx_t idx) const {
    return index_not_found(idx) ? invalid_index : index_switch_to_index_[idx];
  }

  bool operator==(ElectronSymmetric const &rhs) const;
  bool operator!=(ElectronSymmetric const &rhs) const;

  std::map<bit_t, std::pair<idx_t, idx_t>> ups_lower_upper_;
  std::vector<bit_t> dns_;

  std::map<bit_t, std::pair<idx_t, idx_t>> dns_lower_upper_;
  std::vector<bit_t> ups_;

  std::vector<complex> character_switch_;

// private:
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

  std::vector<double> norms_;
  std::vector<idx_t> index_switch_to_index_;

  idx_t size_;
};

} // namespace hydra
