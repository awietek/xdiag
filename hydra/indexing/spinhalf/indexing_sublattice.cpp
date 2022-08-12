#include "indexing_sublattice.h"

#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/group_action/group_action_operations.h>
#include <hydra/symmetries/symmetry_operations.h>
#include <hydra/utils/logger.h>

#include <algorithm>

namespace hydra::indexing::spinhalf {

template <typename bit_t>
void compute_rep_search_range(
    std::vector<bit_t> const &reps, int n_postfix_bits,
    std::unordered_map<bit_t, gsl::span<bit_t const>> &rep_search_range) {
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

  // find all representatives
  for (auto state : combinatorics::Subsets<bit_t>(n_sites)) {
    bit_t rep = representative(state);
    if (state == rep) {
      double norm = symmetries::norm(rep, group_action_, irrep);
      if (std::abs(norm) > 1e-6) {
        reps_.push_back(rep);
        norms_.push_back(norm);
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

  // find all representatives
  for (auto state : combinatorics::Combinations<bit_t>(n_sites, n_up)) {
    bit_t rep = representative(state);
    if (state == rep) {
      double norm = symmetries::norm(rep, group_action_, irrep);
      if (std::abs(norm) > 1e-6) {
        reps_.push_back(rep);
        norms_.push_back(norm);
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
