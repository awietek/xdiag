#include "fermi_sign.h"

#include <hydra/bitops/bitops.h>

namespace hydra::symmetries {

template <class bit_t>
bool fermi_bool_of_permutation(bit_t state, const int *permutation, int *work) {
  int n_fermion = 0;
  for (int site = 0; state; ++site) {
    if (state & 1) {
      work[n_fermion++] = permutation[site];
    }
    state >>= 1;
  }

  int cnt = 0;
  for (int i = 0; i < n_fermion; i++)
    for (int j = i + 1; j < n_fermion; j++)
      if (work[i] > work[j])
        cnt++;

  return (bool)(cnt & 1);
}

template <class bit_t>
double fermi_sign_of_permutation(bit_t state, const int *permutation,
                                 int *work) {
  return fermi_bool_of_permutation(state, permutation, work) ? -1.0 : 1.0;
}

template <class bit_t>
double fermi_sign_of_permutation_sort(
    bit_t state, const int *permutation,
    int *work) { // "work" needs to be allocated of size 2*n_sites
  int *iota = work;
  int *to = work + bitops::popcnt(state);

  // find out where fermions are mapped to
  int n_fermion = 0;
  for (int site = 0; state; ++site) {
    if (state & 1) {
      iota[n_fermion] = n_fermion;
      to[n_fermion++] = permutation[site];
    }
    state >>= 1;
  }

  // Find sorting permutation -> iota
  std::sort(iota, iota + n_fermion,
            [&to](const int &a, const int &b) { return to[a] < to[b]; });

  // compute sign in O(n_fermions) by cycle decomposition
  double sign = 1.0;
  bit_t visited = ((bit_t)1 << n_fermion) - 1;
  int next, L;
  bit_t mask;
  for (int k = 0; k < n_fermion; ++k) {
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
        sign = -sign;
    }
  }
  return sign;
}

template bool fermi_bool_of_permutation<uint16_t>(uint16_t, const int *, int *);
template bool fermi_bool_of_permutation<uint32_t>(uint32_t, const int *, int *);
template bool fermi_bool_of_permutation<uint64_t>(uint64_t, const int *, int *);

template double fermi_sign_of_permutation<uint16_t>(uint16_t, const int *,
                                                    int *);
template double fermi_sign_of_permutation<uint32_t>(uint32_t, const int *,
                                                    int *);
template double fermi_sign_of_permutation<uint64_t>(uint64_t, const int *,
                                                    int *);

template double fermi_sign_of_permutation_sort<uint16_t>(uint16_t, const int *,
                                                         int *);
template double fermi_sign_of_permutation_sort<uint32_t>(uint32_t, const int *,
                                                         int *);
template double fermi_sign_of_permutation_sort<uint64_t>(uint64_t, const int *,
                                                         int *);

} // namespace hydra::utils
