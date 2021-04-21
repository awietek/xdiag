#include "spinflip.h"

namespace hydra {

template <class bit_t>
Spinflip<bit_t>::Spinflip(int n_sites, bool identity_only)
    : n_sites_(n_sites), sitemask_(((bit_t)1 << n_sites) - 1),
      n_symmetries_(identity_only ? 1 : 2) {}

template <class bit_t>
Spinflip<bit_t>
Spinflip<bit_t>::subgroup(std::vector<int> const &symmetry_numbers) const {
  if (symmetry_numbers == std::vector({0, 1}))
    return Spinflip(n_sites);
  else if (symmetry_numbers == std::vector({0}))
    return Spinflip(n_sites, true);
  else
    HydraLog.err("Error creating subgroup of Spinflip: "
                 "invalid symmetry numbers specified");
}

} // namespace hydra
