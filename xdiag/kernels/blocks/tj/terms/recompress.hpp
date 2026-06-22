// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/bits/bitmask.hpp>

namespace xdiag::kernels::tj {

// Compressed-space "slot surgery" for re-mapping the dn string when ups changes.
//
// The dn sector is stored compressed into the non-up sites. When an up operator
// adds/removes an up at a site, that site leaves/joins the dn complement, so the
// dn string gains/loses a slot at that site's compressed rank. Because tJ forbids
// double occupancy, the affected slot is always EMPTY (a site that gains an up
// had no dn; a site that loses an up could not have had a dn there either), so
// these are pure position shifts of the other dn bits -- O(1), portable, and free
// of any parallel bits extract/deposit.

// Remove the slot at compressed rank r; bits above r shift down by one.
template <typename bit_t>
inline bit_t remove_slot(bit_t x, int64_t r, int64_t nsites) {
  bit_t below = bits::bitmask<bit_t>(nsites, r);             // ranks [0, r)
  bit_t above = x & (bits::bitmask<bit_t>(nsites, nsites) ^
                     bits::bitmask<bit_t>(nsites, r + 1));   // ranks [r+1, .)
  return (x & below) | (above >> 1);
}

// Insert an empty slot at compressed rank r; bits at/above r shift up by one.
template <typename bit_t>
inline bit_t insert_zero_slot(bit_t x, int64_t r, int64_t nsites) {
  bit_t below = bits::bitmask<bit_t>(nsites, r);                          // [0,r)
  bit_t above = x & (bits::bitmask<bit_t>(nsites, nsites) ^ below);       // [r,.)
  return (x & below) | (above << 1);
}

} // namespace xdiag::kernels::tj
