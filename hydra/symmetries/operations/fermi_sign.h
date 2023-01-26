#pragma once
#include <hydra/common.h>
#include <hydra/symmetries/permutation.h>

namespace hydra::symmetries {

std::vector<int> fermi_work(int n_sites);

template <class bit_t>
bool fermi_bool_of_permutation(bit_t state, Permutation const &permutation,
                               std::vector<int> &work);
// "work" needs to be allocated of size n_sites + 4

template <class bit_t>
inline bool fermi_bool_of_permutation(bit_t state,
                                      Permutation const &permutation) {
  auto work = fermi_work(permutation.size());
  return fermi_bool_of_permutation(state, permutation, work);
}

template <class bit_t>
double fermi_sign_of_permutation(bit_t state, Permutation const &permutation,
                                 std::vector<int> &work);
// "work" needs to be allocated of size n_sites

template <class bit_t>
inline bool fermi_sign_of_permutation(bit_t state,
                                      Permutation const &permutation) {
  auto work = fermi_work(permutation.size());
  return fermi_sign_of_permutation(state, permutation, work);
}

std::vector<int> fermi_work_sort(int n_sites);

template <class bit_t>
bool fermi_bool_of_permutation_sort(bit_t state, Permutation const &permutation,
                                    std::vector<int> &work);
// "work" needs to be allocated of size 2*n_sites

template <class bit_t>
double fermi_sign_of_permutation_sort(bit_t state,
                                      Permutation const &permutation,
                                      std::vector<int> &work);
// "work" needs to be allocated of size 2*n_sites

} // namespace hydra::symmetries
