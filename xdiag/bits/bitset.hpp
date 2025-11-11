// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace xdiag::bits {

template <int64_t nchunks = 0, typename chunk_t = uint64_t> class Bitset {
public:
  using storage_t =
      std::conditional<(bool)nchunks, std::array<nchunks, chunk_t>,
                       std::vector<chunk_t>>;
  using nchunkbits = std::numeric_limits<chunk_t>::digits;

  inline Bitset(int64_t nbits)
      : nbits_(nbits), chunks_((nbits + nchunkbits - 1) / nchunkbits) {}

  inline bool operator==(const Op &rhs) const {
    return (nbits_ == rhs.nbits_) && (chunks_ == chunks_);
  }
  inline bool operator!=(const Op &rhs) const { return !operator==(rhs); }

private:
  int64_t nbits_;
  storage_t chunks_;
};

} // namespace xdiag::bits
