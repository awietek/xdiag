#pragma once
#include <hydra/common.h>

namespace hydra::symmetries {

std::vector<int> fermi_work(int n_sites);

template <class bit_t>
bool fermi_bool_of_permutation(
    bit_t state, const int *__restrict permutation,
    int *__restrict work); // "work" needs to be allocated of size n_sites + 4

template <class bit_t>
double fermi_sign_of_permutation(
    bit_t state, const int *permutation,
    int *work); // "work" needs to be allocated of size n_sites

std::vector<int> fermi_work_sort(int n_sites);

template <class bit_t>
bool fermi_bool_of_permutation_sort(
    bit_t state, const int *permutation,
    int *work); // "work" needs to be allocated of size 2*n_sites

template <class bit_t>
double fermi_sign_of_permutation_sort(
    bit_t state, const int *permutation,
    int *work); // "work" needs to be allocated of size 2*n_sites

} // namespace hydra::symmetries
