#pragma once

#include <hydra/common.h>

namespace hydra::utils {

template <class bit_t>
double fermi_sign_of_permutation(
    bit_t state, const int *permutation,
    int *work); // "work" needs to be allocated of size n_sites

template <class bit_t>
double fermi_sign_of_permutation_sort(
    bit_t state, const int *permutation,
    int *work); // "work" needs to be allocated of size 2*n_sites

} // namespace hydra::utils
