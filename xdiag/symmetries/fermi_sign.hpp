// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <type_traits>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/symmetries/permutation.hpp>

namespace xdiag::symmetries {

// Backends for which fermi_bool_of_permutation (and hence the fermionic
// representative-table path) is compiled: native integers and Bitset, which
// provide the bitwise operators it uses. The BitArray backends -- used only by
// the bosonic Schaefer / multiset / partition bases, which never build a
// fermionic table -- are excluded, so the fermionic path is not instantiated
// for them.
template <typename T> struct fermi_capable : std::is_integral<T> {};
template <typename chunk_t, int64_t nchunks>
struct fermi_capable<bits::Bitset<chunk_t, nchunks>> : std::true_type {};

// Whether the site permutation `perm` acting on the fermionic occupation
// `state` picks up a minus sign: applying the permutation reorders the creation
// operators of the occupied sites, and the parity of that reordering (the
// number of inversions among the permuted occupied sites) is returned. `true`
// means sign -1, `false` means sign +1.
//
template <typename bit_t>
bool fermi_bool_of_permutation(bit_t const &state, Permutation const &perm);

} // namespace xdiag::symmetries
