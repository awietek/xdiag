// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

// Compact storage for multiple small integers using bit packing.
//
// BitArray stores multiple small integers (1-8 bits each) packed into
// a single native integer or Bitset. It provides fast get/set operations
// with specialized implementations for both native integers and Bitsets.
//
// Template parameters:
//   bit_t: Storage type - uint32_t, uint64_t, or Bitset<uint64_t, X>
//   nbits: Bits per element (1-8). Determines how many elements fit.
//
// Example:
//   BitArray<uint64_t, 3> arr;  // Stores 21 3-bit values (0-7 each)
//   arr.set(0, 5);              // Set first element to 5
//   int64_t val = arr.get(0);   // Get first element (returns 5)
//
// Storage capacity:
//   maximum_size = total_bits / nbits
//   For uint64_t with nbits=3: 64/3 = 21 elements
//   For Bitset<uint64_t, 2> with nbits=3: 128/3 = 42 elements
template <typename bit_tt, int nbitss> class BitArray {
public:
  using bit_t = bit_tt;
  static constexpr int nbits = nbitss;

  // Maximum number of elements that fit in the storage
  static constexpr int64_t maximum_size =
      std::numeric_limits<bit_t>::digits / nbits;

  // Default constructor, initializes all bits to 0
  BitArray() = default;

  // Raw storage access (for use with BitVector)
  explicit BitArray(bit_t raw) noexcept;
  bit_t raw() const noexcept;

  // Get element at index (returns nbits-bit value).
  // Inlined so the compiler sees the compile-time constants nbits and
  // slot_mask and can eliminate the multiply idx*nbits via strength reduction.
  inline int64_t get(int64_t idx) const noexcept {
    if constexpr (std::is_integral_v<bit_t>) {
      constexpr bit_t slot_mask = bitmask<bit_t>(nbits);
      return (bits_ >> (static_cast<int>(idx) * nbits)) & slot_mask;
    } else {
      return bits_.get_range(idx * nbits, nbits);
    }
  }

  // Set element at index to value (stores lower nbits of value).
  // Inlined for the same reason as get(). Masks value before shifting to
  // avoid a second mask after the shift, and uses bit_t arithmetic throughout
  // to stay within the storage type's width.
  inline void set(int64_t idx, int64_t value) noexcept {
    if constexpr (std::is_integral_v<bit_t>) {
      constexpr bit_t slot_mask = bitmask<bit_t>(nbits);
      const int shift = static_cast<int>(idx) * nbits;
      bits_ &= ~(slot_mask << shift);               // clear slot
      bits_ |= (bit_t(value) & slot_mask) << shift; // write value
    } else {
      bits_.set_range(idx * nbits, nbits, (typename bit_t::chunk_t)value);
    }
  }

  bool operator==(BitArray<bit_t, nbits> const &rhs) const noexcept;
  bool operator!=(BitArray<bit_t, nbits> const &rhs) const noexcept;
  bool operator<(BitArray<bit_t, nbits> const &rhs) const noexcept;
  bool operator<=(BitArray<bit_t, nbits> const &rhs) const noexcept;
  bool operator>(BitArray<bit_t, nbits> const &rhs) const noexcept;
  bool operator>=(BitArray<bit_t, nbits> const &rhs) const noexcept;

private:
  bit_t bits_{};
};

// Create a BitArray with proper initialization for n_slots elements.
// For integral/static bit_t: default construct (storage is fixed-size).
// For BitsetDynamic bit_t: allocate n_slots * nbits bits initialized to zero.
template <typename bitarray_t>
inline bitarray_t make_bitarray(int64_t n_slots) {
  using bit_t = typename bitarray_t::bit_t;
  if constexpr (std::is_same_v<bit_t, BitsetDynamic>) {
    return bitarray_t(bit_t(n_slots * bitarray_t::nbits));
  } else {
    return bitarray_t{};
  }
}

namespace detail {
template <typename T> struct is_bitarray : std::false_type {};
template <typename bit_t, int N>
struct is_bitarray<BitArray<bit_t, N>> : std::true_type {};
template <typename T>
inline constexpr bool is_bitarray_v = is_bitarray<T>::value;
} // namespace detail

template <typename bit_t, int nbits>
std::string to_string(BitArray<bit_t, nbits> const &bits,
                      int64_t size = BitArray<bit_t, nbits>::maximum_size,
                      bool reverse = true);

template <typename bit_t, int nbits>
std::ostream &operator<<(std::ostream &out, BitArray<bit_t, nbits> const &bits);

using BitArray1 = BitArray<uint64_t, 1>;
using BitArray2 = BitArray<uint64_t, 2>;
using BitArray3 = BitArray<uint64_t, 3>;
using BitArray4 = BitArray<uint64_t, 4>;
using BitArray5 = BitArray<uint64_t, 5>;
using BitArray6 = BitArray<uint64_t, 6>;
using BitArray7 = BitArray<uint64_t, 7>;
using BitArray8 = BitArray<uint64_t, 8>;

using BitArrayLong1 = BitArray<BitsetDynamic, 1>;
using BitArrayLong2 = BitArray<BitsetDynamic, 2>;
using BitArrayLong3 = BitArray<BitsetDynamic, 3>;
using BitArrayLong4 = BitArray<BitsetDynamic, 4>;
using BitArrayLong5 = BitArray<BitsetDynamic, 5>;
using BitArrayLong6 = BitArray<BitsetDynamic, 6>;
using BitArrayLong7 = BitArray<BitsetDynamic, 7>;
using BitArrayLong8 = BitArray<BitsetDynamic, 8>;

} // namespace xdiag::bits
