#pragma once

#include <vector>
#include <hydra/utils/openmp_utils.h>
#include <hydra/symmetries/operations/fermi_sign.h>

namespace hydra::symmetries {

// Creates the table of Fermi signs for a given particle number
template <class States, class GroupAction>
std::vector<bool> fermi_bool_table(States &&states,
                                   GroupAction &&group_action) {
  assert(states.n() == group_action.n_sites());

  int n_sites = group_action.n_sites();
  int n_symmetries = group_action.n_symmetries();
  idx_t raw_size = states.size();
  std::vector<bool> fermi_bool_table(raw_size * n_symmetries);
  const int *sym_ptr = group_action.permutation_array().data();

#ifndef _OPENMP
  auto fermi_work = symmetries::fermi_work(n_sites);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    idx_t idx = 0;
    for (auto state : states) {
      fermi_bool_table[sym * raw_size + idx] =
          fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
      ++idx;
    }
    sym_ptr += n_sites;
  }
#else

  for (int sym = 0; sym < n_symmetries; ++sym) {

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
        bool fermi_bool =
            fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
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
