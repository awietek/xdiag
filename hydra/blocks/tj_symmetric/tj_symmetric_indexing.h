#pragma once

#include <map>
#include <vector>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/representation.h>
#include <hydra/symmetries/symmetry_operations.h>

#include <hydra/blocks/tj_symmetric/tj_symmetric_matrix.h>

namespace hydra::indexing {

template <typename bit_t> class tJSymmetricIndexing {
public:
  tJSymmetricIndexing(int n_sites, int nup, int ndn,
                      PermutationGroup permutation_group, Representation irrep);

  PermutationGroupLookup<bit_t> const &group_action() const {
    return group_action_;
  }
  Representation const &irrep() const { return irrep_; }
  idx_t size() const { return size_; }

  idx_t n_reps_up() const { return reps_up_.size(); }
  idx_t n_reps_dn() const { return reps_dn_.size(); }

  bit_t rep_up(idx_t idx_up) const { return reps_up_.at(idx_up); }
  bit_t rep_dn(idx_t idx_dn) const { return reps_dn_.at(idx_dn); }

  idx_t up_offset(idx_t idx_up) const { return up_offsets_.at(idx_up); }
  idx_t dn_offset(idx_t idx_dn) const { return up_offsets_.at(idx_dn); }

  bool fermi_bool_up(int sym, bit_t ups) const {
    return fermi_bool_ups_table_[sym * raw_ups_size_ +
                                 lintable_ups_.index(ups)];
  }

  bool fermi_bool_dn(int sym, bit_t dns) const {
    return fermi_bool_dns_table_[sym * raw_dns_size_ +
                                 lintable_dns_.index(dns)];
  }

  std::pair<idx_t, bool> index_dn_fermi(bit_t dns, int sym, bit_t not_ups,
                                        bit_t fermimask) const {
    bit_t dns_rep = group_action_.apply(sym, dns);
    bit_t dns_rep_c = bitops::extract(dns_rep, not_ups);

    int n_sites = 7;
    lila::Log("not_ups: {}", BSTR(not_ups));
    lila::Log("dns_rep: {}", BSTR(dns_rep));
    lila::Log("dns_rep_c: {}", BSTR(dns_rep_c));

    idx_t idx_dns_rep = lintable_dnsc_.index(dns_rep_c);
    bool fermi_dn = (bitops::popcnt(dns & fermimask) & 1);
    fermi_dn ^= fermi_bool_dn(sym, dns);
    return {idx_dns_rep, fermi_dn};
  }

  std::tuple<idx_t, bool, int>
  index_dn_fermi_sym(bit_t dns, std::vector<int> const &syms,
                     std::vector<bit_t> const &dnss_out,
                     bit_t fermimask) const {
    auto [rep_dns, rep_sym] =
        symmetries::representative_sym_subset(dns, group_action_, syms);
    auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
    if ((it != dnss_out.end()) && (*it == rep_dns)) {
      bool fermi_dn = (bitops::popcnt(dns & fermimask) & 1);
      fermi_dn ^= fermi_bool_dn(rep_sym, dns);
      return {std::distance(dnss_out.begin(), it), fermi_dn, rep_sym};
    } else {
      return {invalid_index, false, rep_sym};
    }
  }

  std::pair<idx_t, idx_t> sym_limits_up(bit_t ups) const {
    return sym_limits_up_.at(lintable_ups_.index(ups));
  }
  std::pair<idx_t, idx_t> sym_limits_dn(bit_t dns) const {
    return sym_limits_dn_.at(lintable_dns_.index(dns));
  }

  std::vector<int> syms_up(bit_t ups) const {
    auto [l, u] = sym_limits_up(ups);
    return std::vector<int>(syms_up_.begin() + l, syms_up_.begin() + u);
  }
  std::vector<int> syms_dn(bit_t dns) const {
    auto [l, u] = sym_limits_dn(dns);
    return std::vector<int>(syms_dn_.begin() + l, syms_dn_.begin() + u);
  }

  idx_t index_up(bit_t ups) const {
    return idces_up_[lintable_ups_.index(ups)];
  }

  idx_t index_dn(bit_t dns) const {
    return idces_dn_[lintable_dns_.index(dns)];
  }

  std::vector<bit_t> const &dns_for_up_rep(bit_t ups) const {
    auto [sym_lower, sym_upper] = sym_limits_up(ups);

    // if trivial stabilizer, dns are compressed and need to be deposited
    if ((sym_upper - sym_lower) == 1) {
      return dns_full_;
    } else {
      return dns_for_up_rep_.at(ups);
    }
  }

  std::vector<double> const &norms_for_up_rep(bit_t ups) const {
    auto [sym_lower, sym_upper] = sym_limits_up(ups);
    if ((sym_upper - sym_lower) == 1) {
      return norms_dns_full_;
    } else {
      return norms_for_up_rep_.at(ups);
    }
  }

private:
  int n_sites_;
  PermutationGroupLookup<bit_t> group_action_;
  Representation irrep_;

  idx_t raw_ups_size_;
  idx_t raw_dns_size_;
  idx_t raw_upsc_size_;
  idx_t raw_dnsc_size_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dns_;
  indexing::LinTable<bit_t> lintable_upsc_;
  indexing::LinTable<bit_t> lintable_dnsc_;

  std::vector<bool> fermi_bool_ups_table_;
  std::vector<bool> fermi_bool_dns_table_;

  std::vector<bit_t> reps_up_;
  std::vector<idx_t> idces_up_;
  std::vector<int> syms_up_;
  std::vector<std::pair<idx_t, idx_t>> sym_limits_up_;

  std::vector<bit_t> reps_dn_;
  std::vector<idx_t> idces_dn_;
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

} // namespace hydra::indexing
