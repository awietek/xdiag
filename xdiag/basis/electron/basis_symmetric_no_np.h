#pragma once

#include <map>
#include <vector>

#include <xdiag/extern/gsl/span>

#include <xdiag/common.h>

#include <xdiag/combinatorics/fermi_table.h>
#include <xdiag/combinatorics/subsets_indexing.h>

#include <xdiag/symmetries/group_action/group_action_lookup.h>
#include <xdiag/symmetries/operations/group_action_operations.h>
#include <xdiag/symmetries/operations/symmetry_operations.h>
#include <xdiag/symmetries/permutation_group.h>
#include <xdiag/symmetries/representation.h>

namespace xdiag::basis::electron {

template <typename bit_t> class BasisSymmetricNoNp {
public:
  using span_size_t = gsl::span<int64_t const>::size_type;
  using bit_type = bit_t;

  BasisSymmetricNoNp(int64_t n_sites, PermutationGroup permutation_group,
                     Representation irrep);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;

  int64_t dim() const;
  int64_t size() const;

  GroupActionLookup<bit_t> const &group_action() const;
  Representation const &irrep() const;

  inline int64_t n_rep_ups() const { return reps_up_.size(); }
  inline bit_t rep_ups(int64_t idx_ups) const { return reps_up_[idx_ups]; }
  inline int64_t ups_offset(int64_t idx_ups) const {
    return ups_offset_[idx_ups];
  }

  // index and fermi sign for dns with trivial stabilizer
  inline int64_t index_dns(bit_t dns) const { return lintable_dns_.index(dns); }

  inline std::pair<int64_t, bool> index_dns_fermi(bit_t dns,
                                                  int64_t sym) const {
    bit_t dns_rep = group_action_.apply(sym, dns);
    int64_t idx_dns_rep = lintable_dns_.index(dns_rep);
    bool fermi_dns = fermi_table_.sign(sym, dns);
    return {idx_dns_rep, fermi_dns};
  }

  inline std::pair<int64_t, bool> index_dns_fermi(bit_t dns, int64_t sym,
                                                  bit_t fermimask) const {
    bit_t dns_rep = group_action_.apply(sym, dns);
    int64_t idx_dns_rep = lintable_dns_.index(dns_rep);
    bool fermi_dns = (bits::popcnt(dns & fermimask) & 1);
    fermi_dns ^= fermi_table_.sign(sym, dns);
    return {idx_dns_rep, fermi_dns};
  }

  // index and fermi sign for dns with non-trivial stabilizer
  inline std::tuple<int64_t, bool, int64_t>
  index_dns_fermi_sym(bit_t dns, gsl::span<int64_t const> syms,
                      gsl::span<bit_t const> dnss_out) const {
    auto [rep_dns, rep_sym] =
        symmetries::representative_sym_subset(dns, group_action_, syms);
    auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
    if ((it != dnss_out.end()) && (*it == rep_dns)) {
      bool fermi_dns = fermi_table_.sign(rep_sym, dns);
      return {std::distance(dnss_out.begin(), it), fermi_dns, rep_sym};
    } else {
      return {invalid_index, false, rep_sym};
    }
  }

  inline std::tuple<int64_t, bool, int64_t>
  index_dns_fermi_sym(bit_t dns, gsl::span<int64_t const> syms,
                      gsl::span<bit_t const> dnss_out, bit_t fermimask) const {
    auto [rep_dns, rep_sym] =
        symmetries::representative_sym_subset(dns, group_action_, syms);
    auto it = std::lower_bound(dnss_out.begin(), dnss_out.end(), rep_dns);
    if ((it != dnss_out.end()) && (*it == rep_dns)) {
      bool fermi_dns = (bits::popcnt(dns & fermimask) & 1);
      fermi_dns ^= fermi_table_.sign(rep_sym, dns);
      return {std::distance(dnss_out.begin(), it), fermi_dns, rep_sym};
    } else {
      return {invalid_index, false, rep_sym};
    }
  }

  // Retrieving index of representative and symmetries for ups
  inline int64_t index_ups(bit_t ups) const {
    return idces_up_[lintable_ups_.index(ups)];
  }

  inline gsl::span<int64_t const> syms_ups(bit_t ups) const {
    int64_t idx_ups = lintable_ups_.index(ups);
    auto [start, length] = sym_limits_up_[idx_ups];
    return {syms_up_.data() + start, length};
  }

  inline std::pair<int64_t, gsl::span<int64_t const>>
  index_syms_up(bit_t ups) const {
    int64_t idx_ups = lintable_ups_.index(ups);
    int64_t index = idces_up_[idx_ups];
    auto [start, length] = sym_limits_up_[idx_ups];
    return {index, {syms_up_.data() + start, length}};
  }

  // Retrieving dns states and norms for given up configuration
  inline gsl::span<bit_t const> dns_for_ups_rep(bit_t ups) const {
    int64_t idx_ups = index_ups(ups);
    auto [start, length] = dns_limits_[idx_ups];
    return {dns_storage_.data() + start, length};
  }

  inline gsl::span<double const> norms_for_ups_rep(bit_t ups) const {
    int64_t idx_ups = index_ups(ups);
    auto [start, length] = dns_limits_[idx_ups];
    return {norms_storage_.data() + start, length};
  }

  inline std::pair<gsl::span<bit_t const>, gsl::span<double const>>
  dns_norms_for_up_rep(bit_t ups) const {
    int64_t idx_ups = index_ups(ups);
    auto [start, length] = dns_limits_[idx_ups];
    auto dnss = gsl::span<bit_t const>{dns_storage_.data() + start, length};
    auto norms = gsl::span<double const>{norms_storage_.data() + start, length};
    return {dnss, norms};
  }

  // Fermi sign when applying sym on states
  inline bool fermi_bool_ups(int64_t sym, bit_t ups) const {
    return fermi_table_.sign(sym, ups);
  }
  inline bool fermi_bool_dns(int64_t sym, bit_t dns) const {
    return fermi_table_.sign(sym, dns);
  }

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;
  GroupActionLookup<bit_t> group_action_;
  Representation irrep_;

  int64_t raw_ups_size_;
  int64_t raw_dns_size_;

  combinatorics::SubsetsIndexing<bit_t> lintable_ups_;
  combinatorics::SubsetsIndexing<bit_t> lintable_dns_;
  combinatorics::FermiTableSubsets<bit_t> fermi_table_;

  std::vector<bit_t> reps_up_;
  std::vector<int64_t> idces_up_;
  std::vector<int64_t> syms_up_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_up_;

  std::vector<bit_t> dns_storage_;
  std::vector<double> norms_storage_;
  std::vector<int64_t> ups_offset_;
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits_;

  int64_t size_;
};

} // namespace xdiag::basis::electron
