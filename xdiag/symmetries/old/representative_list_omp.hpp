// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative_list.hpp"

#ifdef _OPENMP

#include <omp.h>

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::symmetries {

using span_size_t = gsl::span<int64_t const>::size_type;

template <typename bit_t, typename T, class StatesIndexing, class GroupAction>
inline std::tuple<
    std::vector<bit_t>, std::vector<int64_t>, std::vector<int64_t>,
    std::vector<std::pair<span_size_t, span_size_t>>, std::vector<double>>
representatives_indices_symmetries_limits_norms_omp(
    StatesIndexing &&states_indexing, GroupAction &&group_action,
    arma::Col<T> const &characters) try {

  int64_t size = states_indexing.size();
  std::vector<int64_t> idces;
  try {
    idces.resize(size, invalid_index);
  } catch (...) {
    XDIAG_THROW("Cannot allocate memory for index array");
  }

  // Compute all representatives
  std::vector<std::vector<bit_t>> reps_thread;
  std::vector<std::vector<double>> norms_thread;

#pragma omp parallel
  {
    int myid = omp_get_thread_num();
    int rank = omp_get_num_threads();

#pragma omp single
    {
      reps_thread.resize(rank);
      norms_thread.resize(rank);
    }

    auto const &states_indexing_thread =
        ThreadStatesIndex(states_indexing.states_indices());


    for (auto [state, idx] : states_indexing_thread) {

      // Alternative variant, but seems slower
    //   auto mybegin = states_indexing.states_indices().begin();
    // auto myend = states_indexing.states_indices().end();
    // for (int i = 0; i < myid; ++i) {
    //   if (mybegin != myend) {
    //     ++mybegin;
    //   }
    // }

    // auto advance = [&](typename decltype(states_indexing.states_indices())::iterator_t &it) {
    //   for (int i = 0; i < rank; ++i) {
    //     if (it != myend) {
    //       ++it;
    //     }
    //   }
    // };
    // for (auto it = mybegin; it != myend; ++it) {
    //   bit_t state = (*it).first;
    //   int64_t idx = (*it).second;



      if (is_representative(state, group_action)) {
        double nrm = symmetries::norm(state, group_action, characters);
        if (std::abs(nrm) > 1e-6) {
          idces[idx] = reps_thread[myid]
                           .size(); // every thread has disjoint set of {idx}
          reps_thread[myid].push_back(state);
          norms_thread[myid].push_back(nrm);
        }
      }
    }
#pragma omp barrier

    int64_t offset = 0;
    for (int id = 0; id < myid; ++id) {
      offset += reps_thread[id].size();
    }

    for (auto [state, idx] : states_indexing_thread) {
      if (idces[idx] != invalid_index) {
        idces[idx] += offset;
      }
    }

  } // pragma omp parallel
  std::vector<bit_t> reps = omp::combine_vectors(reps_thread);
  std::vector<double> norms = omp::combine_vectors(norms_thread);

  // Determine the number of syms yielding the representative for each state
  std::vector<int64_t> n_syms_for_state;
  try {
    n_syms_for_state.resize(size, 0);
  } catch (...) {
    XDIAG_THROW("Cannot allocate memory for n_syms_for_state array");
  }

  int64_t n_reps = reps.size();

#pragma omp parallel for
  for (int64_t rep_idx = 0; rep_idx < n_reps; ++rep_idx) {
    bit_t rep = reps[rep_idx];

    for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
      bit_t state = group_action.apply(sym, rep);
      int64_t idx = states_indexing.index(state);
      idces[idx] = rep_idx; // every thread has disjoint set of {idx}
      ++n_syms_for_state[idx];
    }
  }
  // compute size and allocate syms array
  int64_t n_syms = std::accumulate(n_syms_for_state.begin(),
                                   n_syms_for_state.end(), (int64_t)0);

  std::vector<int64_t> syms;
  try {
    syms.resize(n_syms, 0);
  } catch (...) {
    XDIAG_THROW("Cannot allocate memory for symmetry array");
  }

  // compute the sym offsets
  std::vector<int64_t> n_syms_for_state_offset;
  try {
    n_syms_for_state_offset.resize(size, 0);
  } catch (...) {
    XDIAG_THROW("Cannot allocate memory for n_syms_for_state_offset array");
  }
  std::partial_sum(n_syms_for_state.begin(), n_syms_for_state.end() - 1,
                   n_syms_for_state_offset.begin() + 1);

  // set the sym_limits
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits;
  try {
    sym_limits.resize(size);
  } catch (...) {
    XDIAG_THROW("Cannot allocate memory for symmetry limits array");
  }

#pragma omp parallel for
  for (int64_t idx = 0; idx < size; ++idx) {
    sym_limits[idx] = {n_syms_for_state_offset[idx], n_syms_for_state[idx]};
  }

  // empty n_syms_for_state again, now used as a counter
  std::fill(n_syms_for_state.begin(), n_syms_for_state.end(), 0);

  // calculate the symmetries yielding the representative
  auto const &group = group_action.permutation_group();

#pragma omp parallel for
  for (int64_t rep_idx = 0; rep_idx < n_reps; ++rep_idx) {
    bit_t rep = reps[rep_idx];

    for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
      bit_t state = group_action.apply(sym, rep);
      int64_t idx = states_indexing.index(state);

      int64_t sym_inv = group.inv(sym);
      int64_t idx_sym = n_syms_for_state_offset[idx] + n_syms_for_state[idx]++;
      syms[idx_sym] = sym_inv;
    }
  }

  return {reps, idces, syms, sym_limits, norms};
} XDIAG_CATCH

} // namespace xdiag::symmetries
#endif
