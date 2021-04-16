#include <cassert>
#include <fstream>
#include <iostream>

#include "spacegroup.h"

namespace hydra {

SpaceGroup::SpaceGroup(std::vector<std::vector<int>> const& symmetries)
    : n_sites_(symmetries[0].size()), n_symmetries_(symmetries.size()),
      symmetries_(symmetries), symmetries_internal_(n_sites_ * n_symmetries_) {
  for (int i = 0; i < (int)symmetries.size(); ++i) {
    
    // Check whether lattice symmetries are well-formed
    assert((int)symmetries[i].size() == n_sites_);
    assert(detail::is_valid_permutation(symmetries[i]));
    std::copy(symmetries[i].begin(), symmetries[i].end(),
              symmetries_internal_.begin() + i * n_sites_);
  }
}

SpaceGroup
SpaceGroup::subgroup(std::vector<int> const& symmetry_numbers) const {
  std::vector<std::vector<int>> subgroup_symmetries;
  for (int n_sym : symmetry_numbers) {
    assert((0 <= n_sym) && (n_sym < (int)symmetries_.size()));
    subgroup_symmetries.push_back(symmetries_[n_sym]);
  }
  return SpaceGroup(subgroup_symmetries);
}

SpaceGroup read_spacegroup(std::string filename) {
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

  return SpaceGroup(lattice_symmetries);
}

} // namespace hydra
