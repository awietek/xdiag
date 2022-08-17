#include "fermi_bool_table.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/operations/fermi_sign.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/utils/openmp_utils.h>
#endif

namespace hydra::symmetries {

template <class States>
std::vector<bool> fermi_bool_table_serial(States const &states,
                                          PermutationGroup const &group) {
  int n_sites = group.n_sites();
  int n_symmetries = group.n_symmetries();
  idx_t raw_size = states.size();
  std::vector<bool> fermi_bool_table(raw_size * n_symmetries);
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
  return fermi_bool_table;
}

#ifdef HYDRA_ENABLE_OPENMP
template <class States>
std::vector<bool> fermi_bool_table_omp(States const &states,
                                       PermutationGroup const &group) {
  int n_sites = group.n_sites();
  int n_symmetries = group.n_symmetries();
  idx_t raw_size = states.size();
  std::vector<bool> fermi_bool_table;
  fermi_bool_table.reserve(raw_size * n_symmetries);

  for (int sym = 0; sym < n_symmetries; ++sym) {
    auto const &perm = group[sym];
    std::vector<std::vector<bool>> fermi_bool_table_local;

#pragma omp parallel
    {
      int myid = omp_get_thread_num();
      int rank = omp_get_num_threads();

#pragma omp single
      { fermi_bool_table_local.resize(rank); }
#pragma omp barrier

      auto states_thread = ThreadStates(states);
      // fermi_bool_table_local[myid].resize(states_thread.size());
      auto fermi_work = symmetries::fermi_work(n_sites);
      idx_t idx = 0;
      for (auto state : states_thread) {
        bool fermi_bool = fermi_bool_of_permutation(state, perm, fermi_work);
        fermi_bool_table_local[myid].push_back(fermi_bool);
        ++idx;
      }
#pragma omp barrier

#pragma omp single
      {
        for (int id = 0; id < rank; ++id) {
          fermi_bool_table.insert(fermi_bool_table.end(),
                                  fermi_bool_table_local[id].begin(),
                                  fermi_bool_table_local[id].end());
        }
      }
      //  auto fermi_bool_table_for_sym =
      //       utils::combine_vectors_copy(fermi_bool_table_local);
      //   std::copy(fermi_bool_table_for_sym.begin(),
      //             fermi_bool_table_for_sym.end(),
      //             fermi_bool_table.begin() + sym * raw_size);
      // }
    } // pragma omp parallel
  }

  return fermi_bool_table;
}
#endif

template <class States>
std::vector<bool> fermi_bool_table(States const &states,
                                   PermutationGroup const &group) {
  assert(states.n() == group.n_sites());
#ifdef HYDRA_ENABLE_OPENMP
  return fermi_bool_table_serial(states, group);
#else
  return fermi_bool_table_serial(states, group);
#endif
}

using c16_t = combinatorics::Combinations<uint16_t>;
using c32_t = combinatorics::Combinations<uint32_t>;
using c64_t = combinatorics::Combinations<uint64_t>;

using s16_t = combinatorics::Subsets<uint16_t>;
using s32_t = combinatorics::Subsets<uint32_t>;
using s64_t = combinatorics::Subsets<uint64_t>;

template std::vector<bool> fermi_bool_table<c16_t>(c16_t const &,
                                                   PermutationGroup const &);
template std::vector<bool> fermi_bool_table<c32_t>(c32_t const &,
                                                   PermutationGroup const &);
template std::vector<bool> fermi_bool_table<c64_t>(c64_t const &,
                                                   PermutationGroup const &);

template std::vector<bool> fermi_bool_table<s16_t>(s16_t const &,
                                                   PermutationGroup const &);
template std::vector<bool> fermi_bool_table<s32_t>(s32_t const &,
                                                   PermutationGroup const &);
template std::vector<bool> fermi_bool_table<s64_t>(s64_t const &,
                                                   PermutationGroup const &);

template std::vector<bool>
fermi_bool_table_serial<c16_t>(c16_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_serial<c32_t>(c32_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_serial<c64_t>(c64_t const &, PermutationGroup const &);

template std::vector<bool>
fermi_bool_table_serial<s16_t>(s16_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_serial<s32_t>(s32_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_serial<s64_t>(s64_t const &, PermutationGroup const &);

#ifdef HYDRA_ENABLE_OPENMP
template std::vector<bool>
fermi_bool_table_omp<c16_t>(c16_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_omp<c32_t>(c32_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_omp<c64_t>(c64_t const &, PermutationGroup const &);

template std::vector<bool>
fermi_bool_table_omp<s16_t>(s16_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_omp<s32_t>(s32_t const &, PermutationGroup const &);
template std::vector<bool>
fermi_bool_table_omp<s64_t>(s64_t const &, PermutationGroup const &);
#endif
  
} // namespace hydra::symmetries
