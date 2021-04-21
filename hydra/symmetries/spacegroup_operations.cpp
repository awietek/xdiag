#include "spacegroup_operations.h"

#include <algorithm>
#include <vector>

#include <hydra/utils/bitops.h>

namespace hydra {

std::vector<std::vector<int>> read_permutations(std::string filename) {
  std::vector<std::vector<int>> lattice_symmetries;
  std::ifstream File(filename.c_str());

  if (File.fail()) {
    std::cerr << "Error in read_spacegroup: Could not open file"
              << "with filename [" << filename << "] given. Abort."
              << std::endl;
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
    File >> tobeparsed;
    for (int si = 0; si < n_sites; ++si) {
      int tosite;
      File >> tosite;
      lattice_symmetries[i].push_back(tosite);
    }
  }

  return lattice_symmetries;
}

namespace detail {

bool is_valid_permutation(int n_sites, const int *permutation) {
  for (int i = 0; i < n_sites; ++i) {
    if (std::find(permutation, permutation + n_sites, i) == permutation + n_sites)
      return false;
  }
  return true;
}

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation) {
  bit_t tstate = 0;
  for (int site = 0; site < n_sites; ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16 apply_permutation<uint16>(uint16, int, const int *);
template uint32 apply_permutation<uint32>(uint32, int, const int *);
template uint64 apply_permutation<uint64>(uint64, int, const int *);

template <class bit_t>
double fermi_sign(bit_t state, int n_sites, const int *permutation) {
  std::vector<int> sort(n_sites); // Bad for performance
  int sum = 0;
  int nf, k;

  for (nf = 0, k = 0; k < n_sites; ++k)
    if (utils::gbit(state, k)) {
      sort[nf] = permutation[k];
      ++nf;
    }

  int old_sum;
  while (true) {
    old_sum = sum;
    for (k = 0; k < (nf - 1); ++k)
      if (sort[k + 1] < sort[k]) {
        sum++;
        std::swap(sort[k + 1], sort[k]);
      }
    if (old_sum == sum)
      break;
  }
  return (sum % 2 ? -1 : 1);
}

template double fermi_sign<uint16>(uint16, int, const int *);
template double fermi_sign<uint32>(uint32, int, const int *);
template double fermi_sign<uint64>(uint64, int, const int *);

} // namespace detail
} // namespace hydra
