#include "indexing_sublattice.h"

#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/operations/group_action_operations.h>
#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

#include <algorithm>

namespace hydra::indexing::spinhalf {

template <typename bit_t>
void compute_rep_search_range(
    std::vector<bit_t> const &reps, int n_postfix_bits,
    // std::unordered_map<bit_t, gsl::span<bit_t const>> &rep_search_range) {
    ska::flat_hash_map<bit_t, gsl::span<bit_t const>> &rep_search_range) {

  if (reps.size() > 0) {
    idx_t start = 0;
    idx_t end = 0;
    bit_t previous_prefix = reps[0] >> n_postfix_bits;
    for (auto rep : reps) {
      bit_t prefix = rep >> n_postfix_bits;
      if (prefix != previous_prefix) {
        rep_search_range[previous_prefix] =
            gsl::span<bit_t const>(reps.data() + start, end - start);
        previous_prefix = prefix;
        start = end;
      }
      ++end;
    }
    rep_search_range[previous_prefix] =
        gsl::span<bit_t const>(reps.data() + start, end - start);
  }
}

template <typename bit_t, int n_sublat>
IndexingSublattice<bit_t, n_sublat>::IndexingSublattice(
    int n_sites, PermutationGroup permutation_group, Representation irrep)
    : n_sites_(n_sites), sz_conserved_(false), n_up_(invalid_n),
      n_postfix_bits_(n_sites - std::min(maximum_prefix_bits, n_sites)),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep) {
  using combinatorics::Subsets;

  int n_sites_sublat = n_sites_ / n_sublat;
  int n_leading = n_sites_sublat;
  int n_trailing = (n_sublat - 1) * n_sites_sublat;

  for (auto prefix : Subsets<bit_t>(n_leading)) {

    // if prefix is not rep, the full state also cannot be rep
    bit_t prefix_rep = group_action_.reps_[n_sublat-1][prefix];
    if (prefix_rep < prefix) {
      continue;
    }

    for (auto postfix : Subsets<bit_t>(n_trailing)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, group_action_, irrep);
        if (std::abs(norm) > 1e-6) {
          reps_.push_back(rep);
          norms_.push_back(norm);
        }
      }
    }
  }

  reps_.shrink_to_fit();
  norms_.shrink_to_fit();
  compute_rep_search_range(reps_, n_postfix_bits_, rep_search_range_);
  begin_ = iterator_t(reps_, 0);
  end_ = iterator_t(reps_, reps_.size());
}

template <typename bit_t, int n_sublat>
IndexingSublattice<bit_t, n_sublat>::IndexingSublattice(
    int n_sites, int n_up, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_postfix_bits_(n_sites - std::min(maximum_prefix_bits, n_sites)),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep) {
  using bitops::popcnt;
  using combinatorics::Combinations;
  using combinatorics::Subsets;

  int n_sites_sublat = n_sites_ / n_sublat;
  int n_leading = n_sites_sublat;
  int n_trailing = (n_sublat - 1) * n_sites_sublat;

  for (auto prefix : Subsets<bit_t>(n_leading)) {

    // if prefix is incompatible with n_up, continue
    int n_up_prefix = popcnt(prefix);
    int n_up_postfix = n_up - n_up_prefix;
    if ((n_up_postfix < 0) || (n_up_postfix > n_trailing)) {
      continue;
    }

    // if prefix is not rep, the full state also cannot be rep
    bit_t prefix_rep = group_action_.reps_[n_sublat-1][prefix];
    if (prefix_rep < prefix) {
      continue;      
    }

    for (auto postfix : Combinations<bit_t>(n_trailing, n_up_postfix)) {
      bit_t state = (prefix << n_trailing) | postfix;
      bit_t rep = representative(state);
      if (state == rep) {
        double norm = symmetries::norm(rep, group_action_, irrep);
        if (std::abs(norm) > 1e-6) {
          reps_.push_back(rep);
          norms_.push_back(norm);
        }
      }
    }
  }

  reps_.shrink_to_fit();
  norms_.shrink_to_fit();
  compute_rep_search_range(reps_, n_postfix_bits_, rep_search_range_);
  begin_ = iterator_t(reps_, 0);
  end_ = iterator_t(reps_, reps_.size());
}

template <typename bit_t, int n_sublat>
idx_t IndexingSublattice<bit_t, n_sublat>::index_of_representative(
    bit_t rep) const {
  bit_t prefix = rep >> n_postfix_bits_;
  auto itr = rep_search_range_.find(prefix);
  if (itr == rep_search_range_.end()) {
    return invalid_index;
  } else {
    gsl::span<bit_t const> search_range = itr->second;
    auto it = std::lower_bound(search_range.begin(), search_range.end(), rep);
    return (*it == rep) ? &(*it) - reps_.data() : invalid_index;
  }
}

template <typename bit_t, int n_sublat>
idx_t IndexingSublattice<bit_t, n_sublat>::index(bit_t state) const {
  bit_t rep = representative(state);
  return index_of_representative(rep);
}

template <typename bit_t, int n_sublat>
std::pair<idx_t, int>
IndexingSublattice<bit_t, n_sublat>::index_sym(bit_t state) const {
  auto [rep, sym] = group_action_.representative_sym(state);
  idx_t idx = index_of_representative(rep);
  return {idx, sym};
}

template <typename bit_t, int n_sublat>
std::pair<idx_t, gsl::span<int const>>
IndexingSublattice<bit_t, n_sublat>::index_syms(bit_t state) const {
  auto [rep, syms] = group_action_.representative_syms(state);
  idx_t idx = index_of_representative(rep);
  return {idx, syms};
}

template class IndexingSublattice<uint16_t, 1>;
template class IndexingSublattice<uint32_t, 1>;
template class IndexingSublattice<uint64_t, 1>;

template class IndexingSublattice<uint16_t, 2>;
template class IndexingSublattice<uint32_t, 2>;
template class IndexingSublattice<uint64_t, 2>;

template class IndexingSublattice<uint16_t, 3>;
template class IndexingSublattice<uint32_t, 3>;
template class IndexingSublattice<uint64_t, 3>;

template class IndexingSublattice<uint16_t, 4>;
template class IndexingSublattice<uint32_t, 4>;
template class IndexingSublattice<uint64_t, 4>;

template class IndexingSublattice<uint16_t, 5>;
template class IndexingSublattice<uint32_t, 5>;
template class IndexingSublattice<uint64_t, 5>;

} // namespace hydra::indexing::spinhalf
