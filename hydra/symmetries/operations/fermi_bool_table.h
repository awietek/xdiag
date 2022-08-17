#pragma once
#include <hydra/symmetries/permutation_group.h>
#include <vector>

namespace hydra::symmetries {

template <class States>
std::vector<bool> fermi_bool_table(States const &states,
                                   PermutationGroup const &group);

template <class States>
std::vector<bool> fermi_bool_table_serial(States const &states,
                                          PermutationGroup const &group);

#ifdef HYDRA_ENABLE_OPENMP
template <class States>
std::vector<bool> fermi_bool_table_omp(States const &states,
                                       PermutationGroup const &group);
#endif

} // namespace hydra::symmetries
