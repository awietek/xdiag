#pragma once

#include <map>
#include <utility>
#include <vector>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/representation.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t = std_bit_t,
          class GroupAction = PermutationGroupLookup<bit_t>>
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

  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline GroupAction const &group_action() const { return group_action_; }
  inline Representation const &irrep() const { return irrep_; }

  inline idx_t size() const { return size_; }

  bool operator==(ElectronSymmetric const &rhs) const;
  bool operator!=(ElectronSymmetric const &rhs) const;


  // inline std::pair<bit_t, bit_t> representative(bit_t ups, bit_t dns)
  // const {
  //   idx_t idx_up = index_up(ups);
  //   auto [sym_lower, sym_upper] = sym_limits_up(ups);
  //   bit_t rep_ups = reps_up_[idx_up];

  //   if (sym_upper - sym_lower == 1) { // trivial stabilizer of ups
  //   (likely)
  //     bit_t rep_dns = group_action_.apply(syms_up_[sym_lower], dns);
  //     return {rep_ups, rep_dns};
  //   } else { // with non-trivial up stabilizer, perform binary search
  //     // Determine dns representative
  //     bit_t rep_dns = std::numeric_limits<bit_t>::max();
  //     for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
  //       bit_t tdns = group_action_.apply(syms_up_[sym_idx], dns);
  //       if (tdns < rep_dns)
  //         rep_dns = tdns;
  //     }
  //     return {rep_ups, rep_dns};
  //   }
  // }

  // inline std::pair<idx_t, complex> index_coeff(bit_t ups, bit_t dns)
  // const {
  //   idx_t idx_up = index_up(ups);
  //   bit_t rep_up = reps_up_[idx_up];
  //   idx_t up_offset = up_offsets_[idx_up];
  //   auto [sym_lower, sym_upper] = sym_limits_up(ups);

  //   if (sym_upper - sym_lower == 1) { // trivial stabilizer of ups
  //   (likely)
  //     int sym = syms_up_[sym_lower];
  //     bit_t rep_dns = group_action_.apply(sym, dns);
  //     complex chi = irrep_.character(sym);
  //     complex coeff =
  //         (fermi_bool_up(sym, ups) == fermi_bool_dn(sym, dns)) ? chi :
  //         -chi;
  //     return {up_offset + lintable_dns_.index(rep_dns), coeff};
  //   } else { // with non-trivial up stabilizer, perform binary search

  //     // Determine dns representative
  //     bit_t rep_dns = std::numeric_limits<bit_t>::max();
  //     int rep_sym = 0;
  //     for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
  //       int sym = syms_up_[sym_idx];
  //       bit_t tdns = group_action_.apply(sym, dns);
  //       if (tdns < rep_dns) {
  //         rep_dns = tdns;
  //         rep_sym = sym;
  //       }
  //     }
  //     auto const &dnss = dns_for_up_rep_.at(rep_up);
  //     auto it = std::lower_bound(dnss.begin(), dnss.end(), rep_dns);
  //     complex chi = irrep_.character(rep_sym);
  //     complex coeff =
  //         (fermi_bool_up(rep_sym, ups) == fermi_bool_dn(rep_sym, dns)) ?
  //         chi
  //                                                                      :
  //                                                                      -chi;
  //     if ((it == dnss.end()) || (*it != rep_dns))
  //       return {invalid_index, 0.};
  //     else
  //       return {up_offset + std::distance(dnss.begin(), it), coeff};
  //     ;
  //   }
  // }



  inline idx_t index(bit_t ups, bit_t dns) const {
    idx_t idx_up = index_up(ups);
    bit_t rep_up = reps_up_[idx_up];
    idx_t up_offset = up_offsets_[idx_up];
    auto [sym_lower, sym_upper] = sym_limits_up(ups);

    if (sym_upper - sym_lower == 1) { // trivial stabilizer of ups (likely)
      bit_t rep_dns = group_action_.apply(syms_up_[sym_lower], dns);
      return up_offset + lintable_dns_.index(rep_dns);
    } else { // with non-trivial up stabilizer, perform binary search

      // Determine dns representative
      bit_t rep_dns = std::numeric_limits<bit_t>::max();
      for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
        bit_t tdns = group_action_.apply(syms_up_[sym_idx], dns);
        if (tdns < rep_dns)
          rep_dns = tdns;
      }
      auto const &dnss = dns_for_up_rep_.at(rep_up);
      auto it = std::lower_bound(dnss.begin(), dnss.end(), rep_dns);
      if ((it == dnss.end()) || (*it != rep_dns)) {
        return invalid_index;
      } else {
        return up_offset + std::distance(dnss.begin(), it);
      }
    }
  }


  inline idx_t index_up(bit_t ups) const {
    return idces_up_[lintable_ups_.index(ups)];
  }

  inline idx_t index_dn(bit_t dns) const {
    return idces_dn_[lintable_dns_.index(dns)];
  }


  inline std::pair<idx_t, idx_t> sym_limits_up(bit_t ups) const {
    return sym_limits_up_[lintable_ups_.index(ups)];
  }

  inline std::pair<idx_t, idx_t> sym_limits_dn(bit_t dns) const {
    return sym_limits_dn_[lintable_dns_.index(dns)];
  }

  inline bool fermi_bool_up(int sym, bit_t ups) const {
    return fermi_bool_ups_table_[sym * raw_ups_size_ +
                                 lintable_ups_.index(ups)];
  }
  // inline double fermi_sign_up(int sym, bit_t ups) const {
  //   return fermi_bool_up(sym, ups) ? -1.0 : 1.0;
  // }

  inline bool fermi_bool_dn(int sym, bit_t dns) const {
    return fermi_bool_dns_table_[sym * raw_dns_size_ +
                                 lintable_dns_.index(dns)];
  }
  // inline double fermi_sign_dn(int sym, bit_t dns) const {
  //   return fermi_bool_dn(sym, dns) ? -1.0 : 1.0;
  // }

  inline std::vector<bit_t> const &dns_for_up_rep(bit_t ups) const {
    auto [sym_lower, sym_upper] = sym_limits_up(ups);
    if ((sym_upper - sym_lower) == 1) {
      return dns_full_;
    } else {
      return dns_for_up_rep_.at(ups);
    }
  }

  inline std::vector<double> const &norms_for_up_rep(bit_t ups) const {
    auto [sym_lower, sym_upper] = sym_limits_up(ups);
    if ((sym_upper - sym_lower) == 1) {
      return norms_dns_full_;
    } else {
      return norms_for_up_rep_.at(ups);
    }
  }

  inline std::vector<bit_t> const &ups_for_dn_rep(bit_t dns) const {
    auto [sym_lower, sym_upper] = sym_limits_dn(dns);
    if ((sym_upper - sym_lower) == 1) {
      return ups_full_;
    } else {
      return ups_for_dn_rep_.at(dns);
    }
  }

  inline std::vector<double> const &norms_for_dn_rep(bit_t dns) const {
    auto [sym_lower, sym_upper] = sym_limits_dn(dns);
    if ((sym_upper - sym_lower) == 1) {
      return norms_ups_full_;
    } else {
      return norms_for_dn_rep_.at(dns);
    }
  }

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
