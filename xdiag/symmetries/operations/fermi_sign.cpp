#include "fermi_sign.h"

#include <xdiag/bits/bitops.h>
#include <limits>

namespace xdiag::symmetries {

std::vector<int64_t> fermi_work(int64_t n_sites) {
  return std::vector<int64_t>(n_sites + 4);
}

template <class bit_t>
bool fermi_bool_of_permutation(bit_t state, Permutation const &permutation,
                               std::vector<int64_t> &work) {
  int64_t n_fermion = 0;
  int64_t site = 0;
  while (state) {
    int64_t trailing_zeros = __builtin_ctz(state);
    site += trailing_zeros;
    work[n_fermion++] = permutation[site++];
    state >>= trailing_zeros + 1;
  }

  // fill potentially extraneous work to use loop unrolling
  std::fill(work.data() + n_fermion, work.data() + n_fermion + 4,
            std::numeric_limits<int64_t>::max());

  bool fermi = false;
  for (int64_t i = 0; i < n_fermion; i++) {
    int64_t worki = work[i];

    // Unrolled loop
    for (int64_t j = i + 1; j < n_fermion; j += 4) {
      fermi ^= (worki > work[j]);
      fermi ^= (worki > work[j + 1]);
      fermi ^= (worki > work[j + 2]);
      fermi ^= (worki > work[j + 3]);
    }
  }
  return fermi;

  // Simple implementation:
  // int64_t n_fermion = 0;
  // for (int64_t site = 0; state; ++site) {
  //   if (state & 1) {
  //     work[n_fermion++] = permutation[site];
  //   }
  //   state >>= 1;
  // }

  // int64_t cnt = 0;
  // for (int64_t i = 0; i < n_fermion; i++)
  //   for (int64_t j = i + 1; j < n_fermion; j++)
  //     if (work[i] > work[j])
  //       cnt++;

  // return (bool)(cnt & 1);
}

template <class bit_t>
double fermi_sign_of_permutation(bit_t state, Permutation const &permutation,
                                 std::vector<int64_t> &work) {
  return fermi_bool_of_permutation(state, permutation, work) ? -1.0 : 1.0;
}

std::vector<int64_t> fermi_work_sort(int64_t n_sites) {
  return std::vector<int64_t>(2 * n_sites, 0);
}

template <class bit_t>
bool fermi_bool_of_permutation_sort(bit_t state, Permutation const &permutation,
                                    std::vector<int64_t> &work) {
  // "work" needs to be allocated of size 2*n_sites

  // int64_t n_fermion = 0;
  // int64_t site = 0;
  // // while (state) {
  // //   if (state & 1) {
  // //     work[n_fermion++] = permutation[site];
  // //   }
  // //   state >>= 1;
  // //   ++site;
  // // }

  // // int64_t cnt = 0;
  // // for (int64_t i = 0; i < n_fermion; i++)
  // //   for (int64_t j = i + 1; j < n_fermion; j++)
  // //     if (work[i] > work[j])
  // //       cnt++;

  // return true;

  int64_t *iota = work.data();
  int64_t *to = work.data() + bits::popcnt(state);

  // find out where fermions are mapped to
  int64_t n_fermion = 0;
  for (int64_t site = 0; state; ++site) {
    if (state & 1) {
      iota[n_fermion] = n_fermion;
      to[n_fermion++] = permutation[site];
    }
    state >>= 1;
  }

  // Find sorting permutation -> iota
  std::sort(iota, iota + n_fermion, [&to](const int64_t &a, const int64_t &b) {
    return to[a] < to[b];
  });

  // compute sign in O(n_fermions) by cycle decomposition
  bool sign = false;
  bit_t visited = ((bit_t)1 << n_fermion) - 1;
  int64_t next, L;
  bit_t mask;
  for (int64_t k = 0; k < n_fermion; ++k) {
    if (visited & (1 << k)) {
      next = k;
      L = 0;
      mask = ((bit_t)1 << next);
      while (visited & mask) {
        ++L;
        visited ^= mask;
        next = iota[next];
        mask = ((bit_t)1 << next);
      }
      if (!(L & 1))
        sign = !sign;
    }
  }
  return sign;
}

template <class bit_t>
double fermi_sign_of_permutation_sort(bit_t state,
                                      Permutation const &permutation,
                                      std::vector<int64_t> &work) {
  // "work" needs to be allocated of size 2*n_sites
  return (fermi_bool_of_permutation_sort(state, permutation, work)) ? -1.0
                                                                    : 1.0;
}
template bool fermi_bool_of_permutation<uint16_t>(uint16_t, Permutation const &,
                                                  std::vector<int64_t> &);
template bool fermi_bool_of_permutation<uint32_t>(uint32_t, Permutation const &,
                                                  std::vector<int64_t> &);
template bool fermi_bool_of_permutation<uint64_t>(uint64_t, Permutation const &,
                                                  std::vector<int64_t> &);

template double fermi_sign_of_permutation<uint16_t>(uint16_t,
                                                    Permutation const &,
                                                    std::vector<int64_t> &);
template double fermi_sign_of_permutation<uint32_t>(uint32_t,
                                                    Permutation const &,
                                                    std::vector<int64_t> &);
template double fermi_sign_of_permutation<uint64_t>(uint64_t,
                                                    Permutation const &,
                                                    std::vector<int64_t> &);

template bool fermi_bool_of_permutation_sort<uint16_t>(uint16_t,
                                                       Permutation const &,
                                                       std::vector<int64_t> &);
template bool fermi_bool_of_permutation_sort<uint32_t>(uint32_t,
                                                       Permutation const &,
                                                       std::vector<int64_t> &);
template bool fermi_bool_of_permutation_sort<uint64_t>(uint64_t,
                                                       Permutation const &,
                                                       std::vector<int64_t> &);

template double
fermi_sign_of_permutation_sort<uint16_t>(uint16_t, Permutation const &,
                                         std::vector<int64_t> &);
template double
fermi_sign_of_permutation_sort<uint32_t>(uint32_t, Permutation const &,
                                         std::vector<int64_t> &);
template double
fermi_sign_of_permutation_sort<uint64_t>(uint64_t, Permutation const &,
                                         std::vector<int64_t> &);

} // namespace xdiag::symmetries
