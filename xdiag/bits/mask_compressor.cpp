// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "mask_compressor.hpp"

namespace xdiag::bits {

// One-time precompute of the "compress" move masks, Hacker's Delight 2nd ed.,
// figure 7-10. The inner shift chain is a parallel prefix over the full word
// width (so it is independent of i); the per-round move mask `mv_[i]` is what
// operator() replays cheaply.
template <typename bit_t>
MaskCompressor<bit_t>::MaskCompressor(bit_t mask) : mask_(mask) {
  bit_t m = mask;
  bit_t mk = (bit_t)(~m << 1); // count zeros to the right
  for (int i = 0; i < kIters; ++i) {
    bit_t mp = mk; // parallel prefix over the full word width
    mp ^= (bit_t)(mp << 1);
    mp ^= (bit_t)(mp << 2);
    mp ^= (bit_t)(mp << 4);
    mp ^= (bit_t)(mp << 8);
    mp ^= (bit_t)(mp << 16);
    if constexpr (kBits > 32) {
      mp ^= (bit_t)(mp << 32);
    }
    bit_t mv = mp & m; // bits to move on this round
    mv_[i] = mv;
    m = (bit_t)((m ^ mv) | (mv >> (1 << i)));
    mk = (bit_t)(mk & ~mp);
  }
}

template class MaskCompressor<uint32_t>;
template class MaskCompressor<uint64_t>;

} // namespace xdiag::bits
