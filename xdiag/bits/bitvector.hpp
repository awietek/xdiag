// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>

#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

template <typename bit_t> class BitVectorReference;

template <typename bit_t> class BitVector {
public:
  BitVector() = default;
  BitVector(int64_t nbits, int64_t size);

  inline bit_t operator[](int64_t index) const {
    return storage_.get_bits(index * nbits_, nbits_);
  }
  inline BitVectorReference<bit_t> operator[](int64_t index) {
    return BitVectorReference<bit_t>(nbits_, index, storage_);
  }

  bool operator==(BitVector<bit_t> const &rhs) const;
  bool operator!=(BitVector<bit_t> const &rhs) const;

private:
  int64_t nbits_;
  int64_t size_;
  Bitset<bit_t> storage_; // dynamically sized Bitset
};

template <typename bit_t> class BitVectorReference {
public:
  inline BitVectorReference(int64_t nbits, int64_t index,
                            Bitset<bit_t> &storage)
      : nbits_(nbits), index_(index), storage_(storage) {}
  inline BitVectorReference &operator=(bit_t bits) {
    storage_.set_bits(index_ * nbits_, nbits_, bits);
    return *this;
  }
  inline operator bit_t() const {
    return storage_.get_bits(index_ * nbits_, nbits_);
  }

private:
  int64_t nbits_;
  int64_t index_;
  Bitset<bit_t> &storage_;
};

} // namespace xdiag::bits
