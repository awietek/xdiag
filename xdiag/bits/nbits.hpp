// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/bits/bitarray.hpp>

namespace xdiag::bits {

// Number of bits per site for a given bit storage type.
// Defaults to 1 (plain integers, Bitset); BitArray<bit_t, N> gives N.
template <typename bit_t> constexpr int64_t nbits = 1;

template <typename bit_t, int N>
constexpr int64_t nbits<BitArray<bit_t, N>> = N;

} // namespace xdiag::bits
