#pragma once

namespace hydra::bits {

template <class bit_t> struct half_bit_type_traits;

template <> struct half_bit_type_traits<uint64_t> {
  using type = uint32_t;
};
template <> struct half_bit_type_traits<uint32_t> {
  using type = uint16_t;
};
template <> struct half_bit_type_traits<uint16_t> {
  using type = uint16_t;
};

template <class bit_t>
using half_bit_t = typename half_bit_type_traits<bit_t>::type;

} // namespace hydra::bits
