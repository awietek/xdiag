#pragma once

#include <map>
#include <vector>

#include <extern/gsl/span>

#include <hydra/common.h>
#include <hydra/indexing/fermi_table.h>
#include <hydra/indexing/lin_table.h>

#include <hydra/symmetries/group_action/group_action_lookup.h>
#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::indexing::tj {

template <typename bit_t> class IndexingSymmetricNp {
public:
  using span_size_t = gsl::span<int const>::size_type;

  IndexingSymmetricNp(int n_sites, int nup, int ndn,
                      PermutationGroup permutation_group, Representation irrep);

  int n_sites() const;
  inline int n_up() const;
  inline int n_dn() const;
  GroupActionLookup<bit_t> const &group_action() const;
  Representation const &irrep() const;
  idx_t size() const;

  idx_t n_rep_ups() const;
  bit_t rep_ups(idx_t idx_up) const;
  idx_t ups_offset(idx_t idx_up) const;

  // index and fermi sign for dns with trivial stabilizer
  std::pair<idx_t, bool> index_dns_fermi(bit_t dns, int sym,
                                         bit_t not_ups) const;

  // NEEDED???
  std::pair<idx_t, bool> index_dns_fermi(bit_t dns, int sym, bit_t not_ups,
                                         bit_t fermimask) const;

  // index and fermi sign for dns with non-trivial stabilizer
  std::tuple<idx_t, bool, int>
  index_dns_fermi_sym(bit_t dns, gsl::span<int const> syms,
                      gsl::span<bit_t const> dnss_out) const;

  // NEEDED ??
  std::tuple<idx_t, bool, int>
  index_dns_fermi_sym(bit_t dns, gsl::span<int const> syms,
                      gsl::span<bit_t const> dnss_out, bit_t fermimask) const;

  // Retrieving index of representative and symmetries for ups
  idx_t index_ups(bit_t ups) const;
  gsl::span<int const> syms_ups(bit_t ups) const;
  std::pair<idx_t, gsl::span<int const>> index_syms_up(bit_t ups) const;

  // Retrieving dns states and norms for given up configuration
  gsl::span<bit_t const> dns_for_ups_rep(bit_t ups) const;
  gsl::span<double const> norms_for_ups_rep(bit_t ups) const;
  std::pair<gsl::span<bit_t const>, gsl::span<double const>>
  dns_norms_for_up_rep(bit_t ups) const;

  // Fermi sign when applying sym on states
  bool fermi_bool_ups(int sym, bit_t ups) const;
  bool fermi_bool_dns(int sym, bit_t dns) const;
  idx_t dnsc_index(bit_t dns) const;

private:
  int n_sites_;
  int n_up_;
  int n_dn_;
  GroupActionLookup<bit_t> group_action_;
  Representation irrep_;

  idx_t raw_ups_size_;
  idx_t raw_dns_size_;
  idx_t raw_dnsc_size_;

  indexing::LinTable<bit_t> lintable_ups_;
  indexing::LinTable<bit_t> lintable_dns_;
  indexing::LinTable<bit_t> lintable_dnsc_;

  FermiTableCombinations<bit_t> fermi_table_ups_;
  FermiTableCombinations<bit_t> fermi_table_dns_;

  std::vector<bit_t> reps_up_;
  std::vector<idx_t> idces_up_;
  std::vector<int> syms_up_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_up_;

  std::vector<idx_t> ups_offset_;
  std::vector<std::pair<span_size_t, span_size_t>> dns_limits_;
  std::vector<bit_t> dns_storage_;
  std::vector<double> norms_storage_;

  idx_t size_;
};

} // namespace hydra::indexing::tj
