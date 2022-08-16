#include "sublattice_stability.h"

#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

namespace hydra::symmetries {

bool is_sublattice_permutation(int n_sublat, int sublat,
                               gsl::span<int const> permutation) {
  int n_sites = permutation.size();
  int n_sites_sublat = n_sites / n_sublat;

  // permutation is sublattice stable, if the sublattice
  // (i.e. indices between begin and end) are mapped to
  // most significant indices (n_sites-n_sites_sublat) to n_sites
  int begin = sublat * n_sites_sublat;
  int end = (sublat + 1) * n_sites_sublat;
  for (int j = begin; j < end; ++j) {
    if (permutation[j] < n_sites - n_sites_sublat) {
      return false;
    }
  }
  return true;
}

int which_sublattice_permutation(int n_sublat,
                                 gsl::span<int const> permutation) {
  int sublat = HYDRA_SUBLATTICE_UNSTABLE;
  for (int s = 0; s < n_sublat; ++s) {
    if (is_sublattice_permutation(n_sublat, s, permutation)) {
      sublat = s;
      break;
    }
  }
  return sublat;
}

bool is_sublattice_stable(int n_sublat, PermutationGroup const &group) {

  for (int sym = 0; sym < group.size(); ++sym) {
    auto permutation = group.permutation(sym);
    int sublat = which_sublattice_permutation(n_sublat, permutation);
    if (sublat == HYDRA_SUBLATTICE_UNSTABLE) {
      return false;
    }
  }
  return true;
}

std::vector<int> sublattice_permutations(int n_sublat, int sublat,
                                         PermutationGroup const &group) {
  std::vector<int> syms;

  for (int sym = 0; sym < group.size(); ++sym) {
    auto permutation = group.permutation(sym);
    if (which_sublattice_permutation(n_sublat, permutation) == sublat) {
      syms.push_back(sym);
    }
  }
  return syms;
}

} // namespace hydra::symmetries
