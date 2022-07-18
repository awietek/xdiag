#include "sublattice_stability.h"

#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::symmetries {

bool is_sublattice_stable(int n_sublat, PermutationGroup const &group) {

  int n_symmetries = group.n_symmetries();
  int n_sites = group.n_sites();

  if ((n_sites % n_sublat) != 0) {
    return false;
  }

  int n_sublat_symmetries[NSubLats];
  get_n_sublat_symmetries<NSubLats>(symmetries, &(n_sublat_symmetries[0]));

  if ((uint32)std::accumulate(n_sublat_symmetries,
                              n_sublat_symmetries + NSubLats,
                              0) != n_symmetries)
    return false;

  // Check if symmetries are in valid order for the sublattice description
  uint32 n_sym = 0;
  for (uint32 sublat = 0; sublat < NSubLats; ++sublat)
    for (int i = 0; i < n_sublat_symmetries[sublat]; ++i, ++n_sym)
      for (uint32 j = sublat * n_sites / NSubLats;
           j < (sublat + 1) * n_sites / NSubLats; ++j)
        if (symmetries[n_sym][j] >= n_sites / NSubLats)
          return false;

  return true;
}


get_sublattice_symmetries(int n_sublat, PermutationGroup const &group){
  
}

  
} // namespace hydra::symmetries
