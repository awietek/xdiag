// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/symmetries/permutation.hpp>

namespace xdiag::symmetries {

// Whether the site permutation `perm` acting on the fermionic occupation
// `state` picks up a minus sign: applying the permutation reorders the creation
// operators of the occupied sites, and the parity of that reordering (the
// number of inversions among the permuted occupied sites) is returned. `true`
// means sign -1, `false` means sign +1.
//
// Generic over bit_t (native integers as well as Bitset / BitArray backends),
// which is why it lives here rather than reusing the integer-only routines in
// symmetries/old/fermi_sign.hpp. It is a header-only template (like the bits::
// helpers it calls) because it is instantiated for every bit_t that appears in
// the representative table. It is only ever evaluated while the table is built
// (the apply path reads the precomputed sign), so the small scratch allocation
// is one-time and negligible next to the surrounding table construction.
template <typename bit_t>
bool fermi_bool_of_permutation(bit_t const &state, Permutation const &perm) {
  int64_t nsites = perm.size();
  std::vector<int64_t> work(nsites);
  int64_t n_fermion = 0;
  for (int64_t site = 0; site < nsites; ++site) {
    if (bits::get(state, site)) {
      work[n_fermion++] = perm[site];
    }
  }
  bool fermi = false;
  for (int64_t i = 0; i < n_fermion; ++i) {
    for (int64_t j = i + 1; j < n_fermion; ++j) {
      fermi ^= (work[i] > work[j]);
    }
  }
  return fermi;
}

} // namespace xdiag::symmetries
