// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <xdiag/common.hpp>
#include <xdiag/symmetries/permutation.hpp>

namespace xdiag::symmetries {

std::vector<int64_t> fermi_work(int64_t nsites);

template <class bit_t>
bool fermi_bool_of_permutation(bit_t state, Permutation const &permutation,
                               std::vector<int64_t> &work);
// "work" needs to be allocated of size nsites + 4

template <class bit_t>
inline bool fermi_bool_of_permutation(bit_t state,
                                      Permutation const &permutation) {
  auto work = fermi_work(permutation.size());
  return fermi_bool_of_permutation(state, permutation, work);
}

template <class bit_t>
double fermi_sign_of_permutation(bit_t state, Permutation const &permutation,
                                 std::vector<int64_t> &work);
// "work" needs to be allocated of size nsites

template <class bit_t>
inline bool fermi_sign_of_permutation(bit_t state,
                                      Permutation const &permutation) {
  auto work = fermi_work(permutation.size());
  return fermi_sign_of_permutation(state, permutation, work);
}

std::vector<int64_t> fermi_work_sort(int64_t nsites);

template <class bit_t>
bool fermi_bool_of_permutation_sort(bit_t state, Permutation const &permutation,
                                    std::vector<int64_t> &work);
// "work" needs to be allocated of size 2*nsites

template <class bit_t>
double fermi_sign_of_permutation_sort(bit_t state,
                                      Permutation const &permutation,
                                      std::vector<int64_t> &work);
// "work" needs to be allocated of size 2*nsites

} // namespace xdiag::symmetries
