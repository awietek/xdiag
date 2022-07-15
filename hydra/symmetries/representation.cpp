#include "representation.h"

#include <fstream>

#include <hydra/utils/logger.h>
#include <numeric>

namespace hydra {

Representation::Representation(std::vector<complex> const &characters)
    : characters_(characters), characters_real_(characters.size()),
      allowed_symmetries_(characters.size(), 0) {
  for (int idx = 0; idx < (int)characters.size(); ++idx) {
    characters_real_[idx] = std::real(characters[idx]);
  }
  std::iota(allowed_symmetries_.begin(), allowed_symmetries_.end(), 0);
}

Representation::Representation(std::vector<complex> const &characters,
                               std::vector<int> const &allowed_symmetries)
    : characters_(characters), characters_real_(characters.size()),
      allowed_symmetries_(allowed_symmetries) {

  assert(characters.size() == allowed_symmetries.size());

  for (int idx = 0; idx < (int)characters.size(); ++idx) {
    characters_real_[idx] = std::real(characters[idx]);
  }
}

Representation
Representation::subgroup(std::vector<int> const &symmetry_numbers) const {
  std::vector<complex> sub_characters;
  for (auto sym : symmetry_numbers)
    sub_characters.push_back(characters_[sym]);
  return Representation(sub_characters);
}

bool Representation::operator==(Representation const &rhs) const {
  return lila::close(lila::Vector(characters_), lila::Vector(rhs.characters_));
}

bool Representation::operator!=(Representation const &rhs) const {
  return !operator==(rhs);
}

Representation read_represenation(std::string filename, std::string repname) {
  std::ifstream File(filename.c_str());
  if (File.fail()) {
    std::cerr << "Error in read_charactertable: Could not open "
              << "file with filename [" << filename << "] given. Abort."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<std::string> names;
  std::vector<std::vector<int>> allowed_symmetries_arr;
  std::vector<std::vector<complex>> characters_arr;

  std::string tobeparsed;
  std::string::size_type pos;

  // Jump to Irreps and parse nreps
  File >> tobeparsed;
  while (tobeparsed.find("[Irreps]") == std::string::npos)
    File >> tobeparsed;
  pos = tobeparsed.find('=');
  int nreps;
  if (pos != std::string::npos)
    nreps = atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
  else
    nreps = -1;
  assert(nreps >= 0);

  // Loop over all representations
  bool found = false;
  for (int i = 0; i < nreps; ++i) {
    File >> tobeparsed;
    while (tobeparsed.find("[Representation]") == std::string::npos)
      File >> tobeparsed;
    pos = tobeparsed.find('=');

    // Get name of representation
    std::string name = tobeparsed.substr(pos + 1, std::string::npos);
    std::vector<int> allowed_ops;
    std::vector<complex> characters;

    // parse number of allowed operations
    while (tobeparsed.find("[AllowedOps]") == std::string::npos)
      File >> tobeparsed;
    pos = tobeparsed.find('=');
    int n_allowed_ops;
    if (pos != std::string::npos)
      n_allowed_ops =
          atoi(tobeparsed.substr(pos + 1, std::string::npos).c_str());
    else
      n_allowed_ops = -1;

    // parse allowed operations
    for (int so = 0; so < n_allowed_ops; ++so) {
      int w;
      File >> w;
      if (!File.good()) {
        std::cerr << "Read Error in read_representation (I)" << std::endl;
        exit(EXIT_FAILURE);
      }
      allowed_ops.push_back(w);
    }

    if (!File.good()) {
      std::cerr << "Read Error in read_representation (II)" << std::endl;
      exit(EXIT_FAILURE);
    }

    // parse bloch factors
    for (int so = 0; so < n_allowed_ops; ++so) {
      double re, im;
      File >> re >> im;
      if (!File.good()) {
        std::cerr << "Read Error in read_representation (III)" << std::endl;
        exit(EXIT_FAILURE);
      }
      characters.push_back(re + complex(0, 1) * im);
    }

    auto rep = Representation(characters, allowed_ops);
    if (name == repname) {
      found = true;
      return rep;
    }
  }
  if (!found)
    Log.err("Error reading representations: "
	    "name not found in file");

  return Representation();
}

bool is_complex(Representation const &rep) {
  for (int i = 0; i < rep.size(); ++i) {
    if (!lila::close(lila::imag(rep.character(i)), 0.0))
      return true;
  }
  return false;
}
bool is_real(Representation const &rep) { return !is_complex(rep); }


Representation TrivialRepresentation(idx_t size) {
  return Representation(std::vector<complex>(size, {1.0, 0.0}));
}

} // namespace hydra
