// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "cdagc_string.hpp"

#include <cstdint>
#include <string>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels {

template <typename bit_t>
CdagCString<bit_t>::CdagCString(int64_t nsites, Monomial const &mono,
                                std::string const &cdag_type,
                                std::string const &c_type) try
    : mask_c_(bits::zero<bit_t>(nsites)),
      mask_cdag_only_(bits::zero<bit_t>(nsites)),
      flipmask_(bits::zero<bit_t>(nsites)),
      signmask_c_(bits::zero<bit_t>(nsites)),
      signmask_cdag_(bits::zero<bit_t>(nsites)) {
  bit_t mask_cdag = bits::zero<bit_t>(nsites);
  int64_t prev_cdag = -1;
  int64_t prev_c = -1;
  bool in_c_block = false;
  for (int64_t k = 0; k < mono.size(); ++k) {
    Op const &op = mono[k];
    std::string type = op.type();
    if (type == cdag_type) {
      // creation-left: no Cdag may follow a C, Cdag sites strictly ascending.
      if (in_c_block) {
        XDIAG_THROW("CdagCString: a Cdag follows a C (the string is not in "
                    "creation-left normal order)");
      }
      int64_t s = op[0];
      if (s <= prev_cdag) {
        XDIAG_THROW("CdagCString: the Cdag block is not strictly ascending");
      }
      prev_cdag = s;
      bits::set(mask_cdag, s);
      signmask_cdag_ ^= bits::bitmask<bit_t>(nsites, s);
    } else if (type == c_type) {
      int64_t s = op[0];
      if (s <= prev_c) {
        XDIAG_THROW("CdagCString: the C block is not strictly ascending");
      }
      in_c_block = true;
      prev_c = s;
      bits::set(mask_c_, s);
      signmask_c_ ^= bits::bitmask<bit_t>(nsites, s);
    } else {
      XDIAG_THROW(std::string("CdagCString: unexpected operator type \"") +
                  type + "\" (expected \"" + cdag_type + "\" or \"" + c_type +
                  "\")");
    }
  }
  // Precompute the state-independent derived masks: pure Cdag sites (excluding
  // the number-operator sites Cdag*C) and the occupation-flip mask.
  mask_cdag_only_ = mask_cdag ^ (mask_cdag & mask_c_);
  flipmask_ = mask_c_ ^ mask_cdag;
}
XDIAG_CATCH

// Explicit instantiations for the bit types the fermion and electron Cdag/C
// string kernels use.
template class CdagCString<uint32_t>;
template class CdagCString<uint64_t>;
template class CdagCString<bits::BitsetDynamic>;
template class CdagCString<bits::BitsetStatic2>;

} // namespace xdiag::kernels
