#include "symmetry_operations.hpp"

#include <algorithm>
#include <fstream>

#include <xdiag/common.hpp>

#include <xdiag/utils/logger.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag {

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
  int64_t n_sites;
  if (pos != std::string::npos)
    n_sites = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    n_sites = -1;

  // Jump to SymmetryOps
  File >> tobeparsed;
  while (tobeparsed.find("[SymmetryOps]") == std::string::npos)
    File >> tobeparsed;

  // Read all symmetries
  int64_t n_symmetries;
  pos = tobeparsed.find('=');
  if (pos != std::string::npos)
    n_symmetries = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    n_symmetries = -1;

  lattice_symmetries.resize(n_symmetries);
  for (int64_t i = 0; i < n_symmetries; ++i) {

    std::vector<int64_t> pv(n_sites);
    File >> tobeparsed;
    for (int64_t si = 0; si < n_sites; ++si) {
      int64_t tosite;
      File >> tosite;
      pv[si] = tosite;
    }
    lattice_symmetries[i] = Permutation(pv);
  }
  return lattice_symmetries;
}

} // namespace xdiag

namespace xdiag::symmetries {

bool is_valid_permutation(int64_t n_sites, const int64_t *permutation) {
  for (int64_t i = 0; i < n_sites; ++i) {
    if (std::find(permutation, permutation + n_sites, i) ==
        permutation + n_sites)
      return false;
  }
  return true;
}

//////////////////////////////////////////////////////////
template <typename bit_t>
bit_t apply_permutation(bit_t state, int64_t n_sites, const int64_t *permutation) {
  bit_t tstate = 0;
  for (int64_t site = 0; site < n_sites; ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, int64_t, const int64_t *);
template uint32_t apply_permutation<uint32_t>(uint32_t, int64_t, const int64_t *);
template uint64_t apply_permutation<uint64_t>(uint64_t, int64_t, const int64_t *);

//////////////////////////////////////////////////////////
template <typename bit_t>
bit_t apply_permutation(bit_t state, gsl::span<int64_t const> permutation) {
  bit_t tstate = 0;
  for (std::size_t site = 0; site < permutation.size(); ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16_t apply_permutation<uint16_t>(uint16_t, gsl::span<int64_t const>);
template uint32_t apply_permutation<uint32_t>(uint32_t, gsl::span<int64_t const>);
template uint64_t apply_permutation<uint64_t>(uint64_t, gsl::span<int64_t const>);

} // namespace xdiag::symmetries
