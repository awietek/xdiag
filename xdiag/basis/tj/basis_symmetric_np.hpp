#pragma once

#include <map>
#include <vector>

#include <xdiag/combinatorics/fermi_table.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/gsl/span>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::basis::tj {

template <typename bit_t> class BasisSymmetricNp {
public:
  using span_size_t = gsl::span<int64_t const>::size_type;
  using bit_type = bit_t;

  BasisSymmetricNp(int64_t n_sites, int64_t nup, int64_t ndn,
                   PermutationGroup permutation_group, Representation irrep);

  int64_t n_sites() const;
  int64_t n_up() const;
  int64_t n_dn() const;

  int64_t dim() const;
  int64_t size() const;

  GroupActionLookup<bit_t> const &group_action() const;
  Representation const &irrep() const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;
  GroupActionLookup<bit_t> group_action_;
  Representation irrep_;

  int64_t raw_ups_size_;
  int64_t raw_dns_size_;
  int64_t raw_dnsc_size_;

  combinatorics::LinTable<bit_t> lintable_ups_;
  combinatorics::LinTable<bit_t> lintable_dns_;
  combinatorics::LinTable<bit_t> lintable_dnsc_;
  combinatorics::FermiTableCombinations<bit_t> fermi_table_ups_;
  combinatorics::FermiTableCombinations<bit_t> fermi_table_dns_;

  std::vector<bit_t> reps_up_;
  std::vector<int64_t> idces_up_;
  std::vector<int64_t> syms_up_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_up_;

  std::vector<int64_t> ups_offset_;
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits_;
  std::vector<bit_t> dns_storage_;
  std::vector<double> norms_storage_;

  int64_t size_;

public:
  int64_t n_rep_ups() const;
  bit_t rep_ups(int64_t idx_up) const;
  int64_t ups_offset(int64_t idx_up) const;

  // index and fermi sign for dns with trivial stabilizer
  std::pair<int64_t, bool> index_dns_fermi(bit_t dns, int64_t sym,
                                           bit_t not_ups) const;

  // NEEDED???
  std::pair<int64_t, bool>
  index_dns_fermi(bit_t dns, int64_t sym, bit_t not_ups, bit_t fermimask) const;

  // index and fermi sign for dns with non-trivial stabilizer
  std::tuple<int64_t, bool, int64_t>
  index_dns_fermi_sym(bit_t dns, gsl::span<int64_t const> syms,
                      gsl::span<bit_t const> dnss_out) const;

  // NEEDED ??
  std::tuple<int64_t, bool, int64_t>
  index_dns_fermi_sym(bit_t dns, gsl::span<int64_t const> syms,
                      gsl::span<bit_t const> dnss_out, bit_t fermimask) const;

  // Retrieving index of representative and symmetries for ups
  int64_t index_ups(bit_t ups) const;
  gsl::span<int64_t const> syms_ups(bit_t ups) const;
  std::pair<int64_t, gsl::span<int64_t const>> index_syms_up(bit_t ups) const;

  // Retrieving dns states and norms for given up configuration
  gsl::span<bit_t const> dns_for_ups_rep(bit_t ups) const;
  gsl::span<double const> norms_for_ups_rep(bit_t ups) const;
  std::pair<gsl::span<bit_t const>, gsl::span<double const>>
  dns_norms_for_up_rep(bit_t ups) const;

  // Fermi sign when applying sym on states
  bool fermi_bool_ups(int64_t sym, bit_t ups) const;
  bool fermi_bool_dns(int64_t sym, bit_t dns) const;
  int64_t dnsc_index(bit_t dns) const;
};

} // namespace xdiag::basis::tj
