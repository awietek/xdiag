#pragma once
#include <utility>
#include <hydra/bitops/bitops.h>

namespace hydra::utils {

template <typename bit_t>
constexpr int tj_site_compressed(bit_t spins, int site) {
  return site - bitops::popcnt(spins & ((bit_t)1 << site));
}

template <typename bit_t>
constexpr std::pair<bit_t, bit_t> upc_dn_to_up_dnc(bit_t upc, bit_t dn) {
  bit_t up = bitops::deposit(upc, ~dn);
  bit_t dnc = bitops::extract(dn, ~up);
  return {up, dnc};
}

template <typename bit_t>
constexpr std::pair<bit_t, bit_t> up_dnc_to_upc_dn(bit_t up, bit_t dnc) {
  bit_t dn = bitops::deposit(dnc, ~up);
  bit_t upc = bitops::extract(up, ~dn);
  return {upc, dn};
}

} // namespace hydra::utils
