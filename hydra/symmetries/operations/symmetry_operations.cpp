#include "symmetry_operations.h"

#include <algorithm>
#include <fstream>

#include <hydra/common.h>

#include <hydra/utils/logger.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/combinatorics/combinations_index.h>

#include <hydra/indexing/lin_table.h>

#include <hydra/symmetries/group_action/group_action.h>
#include <hydra/symmetries/group_action/group_action_lookup.h>

namespace hydra {

std::vector<Permutation> read_permutations(std::string filename) {
  std::vector<Permutation> lattice_symmetries;
  std::ifstream File(filename.c_str());

  if (File.fail()) {
    Log.err("Error in read_spacegroup: Could not open file {}", filename);
    exit(EXIT_FAILURE);
  }

  std::string tobeparsed;
  std::string::size_type pos;

  // Jump to Sites and parse n_sites
  File >> tobeparsed;
  while (tobeparsed.find("[Sites]") == std::string::npos)
    File >> tobeparsed;
  pos = tobeparsed.find('=');
  int n_sites;
  if (pos != std::string::npos)
    n_sites = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    n_sites = -1;

  // Jump to SymmetryOps
  File >> tobeparsed;
  while (tobeparsed.find("[SymmetryOps]") == std::string::npos)
    File >> tobeparsed;

  // Read all symmetries
  int n_symmetries;
  pos = tobeparsed.find('=');
  if (pos != std::string::npos)
    n_symmetries = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    n_symmetries = -1;

  lattice_symmetries.resize(n_symmetries);
  for (int i = 0; i < n_symmetries; ++i) {

    std::vector<int> pv(n_sites);
    File >> tobeparsed;
    for (int si = 0; si < n_sites; ++si) {
      int tosite;
      File >> tosite;
      pv[si] = tosite;
    }
    lattice_symmetries[i] = Permutation(pv);
  }
  return lattice_symmetries;
}

} // namespace hydra

namespace hydra::symmetries {

bool is_valid_permutation(int n_sites, const int *permutation) {
  for (int i = 0; i < n_sites; ++i) {
    if (std::find(permutation, permutation + n_sites, i) ==
        permutation + n_sites)
      return false;
  }
  return true;
}

//////////////////////////////////////////////////////////
template <typename bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation) {
  bit_t tstate = 0;
  for (int site = 0; site < n_sites; ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, int, const int *);
template uint32_t apply_permutation<uint32_t>(uint32_t, int, const int *);
template uint64_t apply_permutation<uint64_t>(uint64_t, int, const int *);

//////////////////////////////////////////////////////////
template <typename bit_t>
bit_t apply_permutation(bit_t state, gsl::span<int const> permutation) {
  bit_t tstate = 0;
  for (std::size_t site = 0; site < permutation.size(); ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, gsl::span<int const>);
template uint32_t apply_permutation<uint32_t>(uint32_t, gsl::span<int const>);
template uint64_t apply_permutation<uint64_t>(uint64_t, gsl::span<int const>);

} // namespace hydra::symmetries
