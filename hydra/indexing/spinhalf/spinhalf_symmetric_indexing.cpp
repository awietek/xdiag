#include "spinhalf_symmetric_indexing.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::indexing {

template <class bit_t>
SpinhalfSymmetricIndexing<bit_t>::SpinhalfSymmetricIndexing(
    int n_sites, int n_up, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), n_up_(n_up),
      group_action_((irrep.allowed_symmetries().size() > 0)
                        ? permutation_group.subgroup(irrep.allowed_symmetries())
                        : permutation_group),
      irrep_(irrep), lin_table_(n_sites, n_up),
      index_for_rep_(combinatorics::binomial(n_sites, n_up),
                     (idx_t)invalid_index),
      sym_limits_for_rep_(combinatorics::binomial(n_sites, n_up)) {
  using combinatorics::Combinations;

  utils::check_nup_spinhalf(n_sites, n_up, "SpinhalfSymmetricIndexing");
  utils::check_n_sites(n_sites, permutation_group);

  // Go through non symmetrized states and register representatives
  idx_t n_reps = 0;
  idx_t idx = 0;
  for (bit_t state : Combinations(n_sites, n_up)) {

    bit_t rep = symmetries::representative(state, group_action_);

    // register state if it's a representative and has non-zero norm
    if (rep == state) {
      double norm = symmetries::norm(rep, group_action_, irrep);
      if (norm > 1e-6) {
        idx_t idx = lin_table_.index(rep);
        reps_.push_back(rep);
        norms_.push_back(norm);
        index_for_rep_[idx] = n_reps++;
      }
    }
    ++idx;
  }

  // Go through non symmetrized states and fill indices for all states
  idx = 0;
  for (bit_t state : Combinations(n_sites, n_up)) {

    bit_t rep = symmetries::representative(state, group_action_);
    idx_t rep_idx = index_for_rep_[lin_table_.index(rep)];
    if (rep_idx != invalid_index) {
      index_for_rep_[idx] = rep_idx;

      // Add syms yielding the representative
      std::vector<int> rep_syms =
          symmetries::mapping_syms(state, rep, group_action_);
      span_size_t start = syms_.size();
      syms_.insert(syms_.end(), rep_syms.begin(), rep_syms.end());
      span_size_t end = syms_.size();
      sym_limits_for_rep_[idx] = {start, end - start};
    }

    ++idx;
  }

  size_ = (idx_t)reps_.size();
}

template class SpinhalfSymmetricIndexing<uint16_t>;
template class SpinhalfSymmetricIndexing<uint32_t>;
template class SpinhalfSymmetricIndexing<uint64_t>;

} // namespace hydra::indexing
