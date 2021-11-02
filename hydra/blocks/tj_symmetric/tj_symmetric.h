#pragma once

#include <map>

#include <hydra/common.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

template <class bit_t = std_bit_t, class GroupAction = PermutationGroupAction>
class tJSymmetric {
public:
  tJSymmetric() = default;
  tJSymmetric(int n_sites, int charge, int sz,
              PermutationGroup permutation_group, Representation irrep);

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

  bool operator==(tJSymmetric const &rhs) const;
  bool operator!=(tJSymmetric const &rhs) const;

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
  int n_symmetries_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dns_;

  idx_t raw_ups_size_;
  idx_t raw_dns_size_;
  std::vector<bool> fermi_bool_ups_table_;
  std::vector<bool> fermi_bool_dns_table_;

  std::vector<idx_t> idces_up_;
  std::vector<bit_t> reps_up_;
  std::vector<int> syms_up_;
  std::vector<std::pair<idx_t, idx_t>> sym_limits_up_;

  std::vector<idx_t> idces_dn_;
  std::vector<bit_t> reps_dn_;
  std::vector<int> syms_dn_;
  std::vector<std::pair<idx_t, idx_t>> sym_limits_dn_;

  std::vector<idx_t> up_offsets_;
  std::vector<bit_t> dns_full_;
  std::vector<double> norms_dns_full_;
  std::map<bit_t, std::vector<bit_t>> dns_for_up_rep_;
  std::map<bit_t, std::vector<double>> norms_for_up_rep_;

  std::vector<idx_t> dn_offsets_;
  std::vector<bit_t> ups_full_;
  std::vector<double> norms_ups_full_;
  std::map<bit_t, std::vector<bit_t>> ups_for_dn_rep_;
  std::map<bit_t, std::vector<double>> norms_for_dn_rep_;

  idx_t size_;
};

} // namespace hydra
