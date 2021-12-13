#pragma once

#include <utility>
#include <vector>

#include <lila/external/gsl/span>

#include <hydra/common.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::symmetries {

using span_size_t = gsl::span<int const>::size_type;

template <typename bit_t, class GroupAction, class LinTable>
std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
           std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits(int npar,
                                          GroupAction const &group_action,
                                          LinTable const &lintable);

template <typename bit_t, class GroupAction, class StateIdxng>
std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
           std::vector<std::pair<span_size_t, span_size_t>>,
           std::vector<double>>
representatives_indices_symmetries_limits_norms(GroupAction &&group_action,
                                                StateIdxng &&idxng,
                                                Representation const &irrep) {
  std::vector<bit_t> reps;
  std::vector<idx_t> idces(idxng.size(), invalid_index);
  std::vector<int> syms;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits(idxng.size());
  std::vector<double> norms;

  // Compute all representatives
  for (auto [state, idx] : idxng.states_indices()) {

    bit_t rep = representative(state, group_action);

    // register state if it's a representative and has non-zero norm
    if (rep == state) {
      double norm = symmetries::norm(rep, group_action, irrep);
      if (norm > 1e-6) {
        idces[idx] = reps.size();
        reps.push_back(rep);
        norms.push_back(norm);
      }
    }

    // Compute indices of up-representatives and stabilizer symmetries
    for (auto [state, idx] : idxng.states_indices()) {

      bit_t rep = representative(state, group_action);
      idx_t rep_idx = idces[idxng.index(rep)];
      if (rep_idx != invalid_index) { // can be invalid if zero norm
        idces[idx] = rep_idx;

        // Add syms yielding the representative
        std::vector<int> rep_syms = mapping_syms(state, rep, group_action);
        span_size_t start = syms.size();
        syms.insert(syms.end(), rep_syms.begin(), rep_syms.end());
        span_size_t end = syms.size();
        sym_limits[idx] = {start, end - start};
      }
    }
  }

  return {reps, idces, syms, sym_limits, norms};
}

} // namespace hydra::symmetries
