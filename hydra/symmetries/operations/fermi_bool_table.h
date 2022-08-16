#pragma once

#include <hydra/symmetries/operations/fermi_sign.h>
#include <hydra/utils/openmp_utils.h>
#include <vector>

namespace hydra::symmetries {

// Creates the table of Fermi signs for a given particle number
template <class States, class GroupAction>
std::vector<bool> fermi_bool_table(States &&states,
                                   GroupAction &&group_action) {
  assert(states.n() == group_action.n_sites());

  int n_sites = group_action.n_sites();
  int n_symmetries = group_action.n_symmetries();
  auto const& group = group_action.permutation_group();
  idx_t raw_size = states.size();
  std::vector<bool> fermi_bool_table(raw_size * n_symmetries);

#ifndef _OPENMP
  auto fermi_work = symmetries::fermi_work(n_sites);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    idx_t idx = 0;
    auto const &perm = group[sym];

    for (auto state : states) {
      fermi_bool_table[sym * raw_size + idx] =
          fermi_bool_of_permutation(state, perm, fermi_work);
      ++idx;
    }
  }
#else

  for (int sym = 0; sym < n_symmetries; ++sym) {
    auto const &perm = group[sym];
    std::vector<std::vector<bool>> fermi_bool_table_local;

    // auto t1 = lila::rightnow();
#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int rank = omp_get_num_threads();

#pragma omp single
      { fermi_bool_table_local.resize(rank); }
#pragma omp barrier

      auto states_thread = ThreadStates(states);
      fermi_bool_table_local[myid].resize(states_thread.size());
      auto fermi_work = symmetries::fermi_work(n_sites);
      idx_t idx = 0;
      for (auto state : states_thread) {
        bool fermi_bool = fermi_bool_of_permutation(state, perm, fermi_work);
        fermi_bool_table_local[myid][idx] = fermi_bool;
        ++idx;
      }
    } // pragma omp parallel
    // lila::timing(t1, lila::rightnow(), "fill");

    // auto t2 = lila::rightnow();
    auto fermi_bool_table_for_sym =
        utils::combine_vectors_copy(fermi_bool_table_local);
    std::copy(fermi_bool_table_for_sym.begin(), fermi_bool_table_for_sym.end(),
              fermi_bool_table.begin() + sym * raw_size);
    // lila::timing(t2, lila::rightnow(), "combine");
    sym_ptr += n_sites;
  }

#endif // ifndef _OPENMP

  return fermi_bool_table;
}

} // namespace hydra::symmetries
