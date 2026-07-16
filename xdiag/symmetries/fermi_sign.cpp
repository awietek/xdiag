// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermi_sign.hpp"

#include <cstdint>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>

namespace xdiag::symmetries {

template <typename bit_t>
bool fermi_bool_of_permutation(bit_t const &state, int64_t const *perm,
                               int64_t nsites) {

  // We walk the sites in increasing order; the inversions contributed by a
  // newly occupied site's image v are exactly the already-seen images greater
  // than v. Tracking the seen images in a bitmask makes that an O(1) popcount
  // per site, so the whole routine is O(nsites) instead of O(nfermion^2).
  bit_t seen = bits::zero<bit_t>(nsites);
  int64_t n_seen = 0;
  bool fermi = false;
  for (int64_t site = 0; site < nsites; ++site) {
    if (bits::get(state, site)) {
      int64_t v = perm[site];
      // already-seen images <= v; those > v are the new inversions
      int64_t n_le = bits::popcount(seen & bits::bitmask<bit_t>(nsites, v + 1));
      fermi ^= (bool)((n_seen - n_le) & 1);
      bits::set(seen, v);
      ++n_seen;
    }
  }
  return fermi;

  // Reference O(nfermion^2) implementation (kept for clarity / cross-checking):
  //   std::vector<int64_t> work(nsites);
  //   int64_t n_fermion = 0;
  //   for (int64_t site = 0; site < nsites; ++site) {
  //     if (bits::get(state, site)) {
  //       work[n_fermion++] = perm[site];
  //     }
  //   }
  //   bool fermi = false;
  //   for (int64_t i = 0; i < n_fermion; ++i) {
  //     for (int64_t j = i + 1; j < n_fermion; ++j) {
  //       fermi ^= (work[i] > work[j]);
  //     }
  //   }
  //   return fermi;
}

template <typename bit_t>
bool fermi_bool_of_permutation(bit_t const &state, Permutation const &perm) {
  return fermi_bool_of_permutation(state, perm.array().data(), perm.size());
}

#define INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(BIT_TYPE)                        \
  template bool fermi_bool_of_permutation(BIT_TYPE const &state,               \
                                          int64_t const *perm,                 \
                                          int64_t nsites);                     \
  template bool fermi_bool_of_permutation(BIT_TYPE const &state,               \
                                          Permutation const &perm);

using namespace bits;
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(uint32_t);
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(uint64_t);
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(BitsetDynamic);
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(BitsetStatic1);
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(BitsetStatic2);
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(BitsetStatic4);
INSTANTIATE_FERMI_BOOL_OF_PERMUTATION(BitsetStatic8);

#undef INSTANTIATE_FERMI_BOOL_OF_PERMUTATION

} // namespace xdiag::symmetries
