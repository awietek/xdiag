#include "representative_list.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/combinations_index.h>

#include <hydra/indexing/lin_table.h>

#include <hydra/symmetries/group_action/group_action.h>
#include <hydra/symmetries/group_action/group_action_lookup.h>
#include <hydra/symmetries/operations/symmetry_operations.h>

#include <omp.h>

namespace hydra::symmetries {

template <typename bit_t, class StatesIndexing, class GroupAction>
std::tuple<std::vector<bit_t>, std::vector<idx_t>, std::vector<int>,
           std::vector<std::pair<span_size_t, span_size_t>>>
representatives_indices_symmetries_limits_norms_omp(
    StatesIndexing &&states_indexing, GroupAction &&group_action,
    Representation const &irrep) {

  idx_t size = states_indexing.size();
  std::vector<idx_t> idces(size, invalid_index);

  // Compute all representatives
  std::vector<std::vector<bit_t>> reps_thread;
  std::vector<std::vector<int>> syms_thread;

#pragma omp parallel
  {
    int myid = omp_get_thread_num();
    int rank = omp_get_num_threads();

    auto states_indexing_thread = ThreadStatesIndex(states_indexing);

#pragma omp single
    { reps_thread.resize(rank); }

    for (auto [state, idx] : states_indexing_thread) {
      bit_t rep = representative(state, group_action);
      if (rep == state) {
        idces[idx] = reps_thread[myid].size(); // just local index
        reps_thread[myid].push_back(rep);
      }
    }

#pragma omp barrier

    idx_t offset = 0;
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

  // Compute indices of up-representatives and stabilizer symmetries
#pragma omp parallel
  {
    int myid = omp_get_thread_num();
    int rank = omp_get_num_threads();

    auto states_indexing_thread = ThreadStatesIndex(states_indexing);

#pragma omp single
    { syms_thread.resize(rank); }

    for (auto [state, idx] : states_indexing_thread) {
      bit_t rep = representative(state, group_action);
      assert(rep <= state);

      idces[idx] = idces[lintable.index(rep)];

      assert(idces[idx] != invalid_index);

      // Determine the symmetries that yield the up-representative
      span_size_t begin = syms_thread[myid].size();
      for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
        if (group_action.apply(sym, state) == rep)
          syms_thread[myid].push_back(sym);
      }
      span_size_t end = syms_thread[myid].size();
      sym_limits[idx] = {begin, end - begin};
    }

#pragma omp barrier

    span_size_t offset = 0;
    for (int id = 0; id < myid; ++id) {
      offset += syms_thread[id].size();
    }

    for (auto [state, idx] : states_indexing_thread) {
      auto [begin, length] = sym_limits[idx];
      sym_limits[idx] = {begin + offset, length};
    }

  } // pragma omp parallel

  // Combine syms
  std::vector<int> syms = utils::combine_vectors(syms_thread);
  std::vector<double> norms;
  return {reps, idces, syms, sym_limits, norms};
}

} // namespace hydra::symmetries
