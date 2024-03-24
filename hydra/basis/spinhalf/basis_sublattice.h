#pragma once

#include <extern/flat_hash_map.h>
#include <extern/gsl/span>

#include <hydra/common.h>

#include <hydra/symmetries/group_action/group_action_sublattice.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::basis::spinhalf {

constexpr int64_t maximum_prefix_bits = 24;

template <typename bit_t, int n_sublat> class BasisSublattice {
  // Public interface
public:
  using iterator_t = typename std::vector<bit_t>::const_iterator;
  using bit_type = bit_t;
  
  BasisSublattice() = default;
  BasisSublattice(int64_t n_sites, PermutationGroup permutation_group,
                  Representation irrep);
  BasisSublattice(int64_t n_sites, int64_t n_up,
                  PermutationGroup permutation_group, Representation irrep);

  iterator_t begin() const;
  iterator_t end() const;
  int64_t size() const;
  int64_t dim() const;

  int64_t index(bit_t state) const;
  inline bit_t state(int64_t idx) const { return reps_[idx]; }
  inline double norm(int64_t idx) const { return norms_[idx]; }

  int64_t n_sites() const;
  bool sz_conserved() const;
  int64_t n_up() const;
  GroupActionSublattice<bit_t, n_sublat> const &group_action() const;
  Representation const &irrep() const;

  bool operator==(BasisSublattice<bit_t, n_sublat> const &rhs) const;
  bool operator!=(BasisSublattice<bit_t, n_sublat> const &rhs) const;

private:
  int64_t n_sites_;
  bool sz_conserved_;
  int64_t n_up_;
  int64_t n_postfix_bits_;

  GroupActionSublattice<bit_t, n_sublat> group_action_;
  Representation irrep_;

  std::vector<bit_t> reps_;
  std::vector<double> norms_;
  // std::unordered_map<bit_t, gsl::span<bit_t const>> rep_search_range_;
  ska::flat_hash_map<bit_t, gsl::span<bit_t const>> rep_search_range_;

  int64_t index_of_representative(bit_t rep) const;

  // functions used in implementation of terms
public:
  inline bit_t representative(bit_t state) const {
    return group_action_.representative(state);
  }
  std::pair<int64_t, int64_t> index_sym(bit_t raw_state) const;
  std::pair<int64_t, gsl::span<int64_t const>>
  index_syms(bit_t raw_state) const;
};

} // namespace hydra::basis::spinhalf
