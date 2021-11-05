#pragma once
#include <hydra/bitops/bitops.h>
#include <utility>

namespace hydra::utils {

template <typename bit_t>
constexpr int tj_site_compressed(bit_t holes, int site) {
  return site - bitops::popcnt(holes & (((bit_t)1 << site) - 1));
}

template <typename bit_t>
constexpr bool has_double_occupations(bit_t up, bit_t dn) {
  return (up & dn);
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

template <typename bit_t>
constexpr std::pair<bit_t, bit_t> holes_spins_to_up_dn(bit_t holes, bit_t spins,
                                                       int n_sites) {
  bit_t mask = ((bit_t)1 << n_sites) - 1;
  bit_t up = bitops::deposit(spins, (bit_t)~holes);
  bit_t dn = bitops::deposit((bit_t)~spins, (bit_t)~holes) & mask;
  return {up, dn};
}

template <typename bit_t>
constexpr std::pair<bit_t, bit_t> up_dn_to_holes_spins(bit_t up, bit_t dn,
                                                       int n_sites) {
  bit_t mask = ((bit_t)1 << n_sites) - 1;
  bit_t holes = (~(up | dn)) & mask;
  bit_t spins = bitops::extract(up, (bit_t)~holes);
  return {holes, spins};
}

} // namespace hydra::utils
