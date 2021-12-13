#include "fermi_bool_table.h"

#include <hydra/combinatorics/combinations_index.h>

#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hydra::symmetries {

template <typename bit_t, class GroupAction>
std::vector<bool> fermi_bool_table(int npar, GroupAction const &group_action) {
  using combinatorics::binomial;
  using combinatorics::CombinationsIndex;
  using combinatorics::CombinationsIndexThread;
  using combinatorics::get_next_pattern;
  using combinatorics::get_nth_pattern;

  int n_sites = group_action.n_sites();
  int n_symmetries = group_action.n_symmetries();
  idx_t raw_size = binomial(n_sites, npar);
  std::vector<bool> fermi_bool_table(raw_size * n_symmetries);
  const int *sym_ptr = group_action.permutation_array().data();

#ifndef _OPENMP
  auto fermi_work = symmetries::fermi_work(n_sites);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    for (auto [state, idx] : CombinationsIndex<bit_t>(n_sites, npar)) {
      fermi_bool_table[sym * raw_size + idx] =
          fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
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

      auto fermi_work = fermi_work(n_sites);

      auto comb = CombinationsIndexThread<bit_t>(n_sites, npar);
      auto begin = comb.begin();
      auto end = comb.end();
      idx_t length = end.idx() - begin.idx();

#pragma omp single
      { fermi_bool_table_local.resize(rank); }

#pragma omp barrier

      fermi_bool_table_local[myid].resize(length);
      auto [state, idxl] = *begin;

      for (idx_t idx = 0; idx < length; ++idx) {
        bool fermi_bool =
            fermi_bool_of_permutation(state, sym_ptr, fermi_work.data());
        fermi_bool_table_local[myid][idx] = fermi_bool;
        state = get_next_pattern(state);
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

template std::vector<bool> fermi_bool_table<uint16_t, PermutationGroupAction>(
    int, PermutationGroupAction const &);
template std::vector<bool> fermi_bool_table<uint32_t, PermutationGroupAction>(
    int, PermutationGroupAction const &);
template std::vector<bool> fermi_bool_table<uint64_t, PermutationGroupAction>(
    int, PermutationGroupAction const &);
template std::vector<bool>
fermi_bool_table<uint16_t, PermutationGroupLookup<uint16_t>>(
    int, PermutationGroupLookup<uint16_t> const &);
template std::vector<bool>
fermi_bool_table<uint32_t, PermutationGroupLookup<uint32_t>>(
    int, PermutationGroupLookup<uint32_t> const &);
template std::vector<bool>
fermi_bool_table<uint64_t, PermutationGroupLookup<uint64_t>>(
    int, PermutationGroupLookup<uint64_t> const &);

} // namespace hydra::symmetries
