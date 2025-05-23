// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::bits {

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

} // namespace xdiag::bits
