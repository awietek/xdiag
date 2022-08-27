#include "fermi_table.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/operations/fermi_sign.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::indexing {

template <class States>
std::vector<bool> init_fermi_table_serial(States const &states,
                                          PermutationGroup const &group) {
  int n_sites = group.n_sites();
  int n_symmetries = group.n_symmetries();
  idx_t raw_size = states.size();
  std::vector<bool> fermi_table(raw_size * n_symmetries);
  auto fermi_work = symmetries::fermi_work(n_sites);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    auto const &perm = group[sym];

    idx_t idx = 0;
    for (auto state : states) {
      fermi_table[sym * raw_size + idx] =
          symmetries::fermi_bool_of_permutation(state, perm, fermi_work);
      ++idx;
    }
  }
  return fermi_table;
}

#ifdef HYDRA_ENABLE_OPENMP
template <class States>
std::vector<bool> init_fermi_table_omp(States const &states,
                                       PermutationGroup const &group) {
  int n_sites = group.n_sites();
  int n_symmetries = group.n_symmetries();
  idx_t raw_size = states.size();
  std::vector<bool> fermi_bool_table;
  fermi_bool_table.reserve(raw_size * n_symmetries);

  for (int sym = 0; sym < n_symmetries; ++sym) {
    auto const &perm = group[sym];
    std::vector<std::vector<bool>> fermi_bool_table_local;

    // Comment: even though every entry in fermi_bool_table would be
    // written by only one thread, std::vector<bool> is special and requires
    // seperate vectors to be built
    
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
        bool fermi_bool =
            symmetries::fermi_bool_of_permutation(state, perm, fermi_work);
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

    } // pragma omp parallel
  }

  return fermi_bool_table;
}
#endif

template <class States>
std::vector<bool> init_fermi_table(States const &states,
                                   PermutationGroup const &group) {
  assert(states.n() == group.n_sites());
#ifdef HYDRA_ENABLE_OPENMP
  return init_fermi_table_omp(states, group);
#else
  return init_fermi_table_serial(states, group);
#endif
}

template <typename bit_t>
FermiTableSubsets<bit_t>::FermiTableSubsets(int n_sites,
                                            PermutationGroup const &group)
    : n_sites_(n_sites),
      table_(init_fermi_table(combinatorics::Subsets<bit_t>(n_sites), group)) {}

template <typename bit_t>
bool FermiTableSubsets<bit_t>::operator==(
    FermiTableSubsets<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (table_ == rhs.table_);
}

template <typename bit_t>
bool FermiTableSubsets<bit_t>::operator!=(
    FermiTableSubsets<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class FermiTableSubsets<uint16_t>;
template class FermiTableSubsets<uint32_t>;
template class FermiTableSubsets<uint64_t>;

template <typename bit_t>
FermiTableCombinations<bit_t>::FermiTableCombinations(
    int n_sites, int n_par, PermutationGroup const &group)
    : raw_size_(combinatorics::binomial(n_sites, n_par)),
      lin_table_(n_sites, n_par),
      table_(init_fermi_table(
          combinatorics::Combinations<bit_t>(n_sites, n_par), group)) {}

template <typename bit_t>
bool FermiTableCombinations<bit_t>::operator==(
    FermiTableCombinations<bit_t> const &rhs) const {
  return (raw_size_ == rhs.raw_size_) && (lin_table_ == rhs.lin_table_) &&
         (table_ == rhs.table_);
}
template <typename bit_t>
bool FermiTableCombinations<bit_t>::operator!=(
    FermiTableCombinations<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class FermiTableCombinations<uint16_t>;
template class FermiTableCombinations<uint32_t>;
template class FermiTableCombinations<uint64_t>;

} // namespace hydra::indexing
