// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cassert>
#include <fstream>
#include <iostream>

#include "charactertable.h"

namespace hydra {

CharacterTable::CharacterTable(
    const SpaceGroup &space_group, const std::vector<std::string> &names,
    const std::vector<std::vector<int>> &allowed_symmetries,
    const std::vector<std::vector<complex>> &characters)
    : space_group_(space_group), names_(names),
      n_symmetries_total_((int)space_group_.n_symmetries()) {
  assert(names.size() == allowed_symmetries.size());
  assert(names.size() == characters.size());

  for (int idx = 0; idx < (int)names.size(); ++idx) {
    std::string name = names[idx];
    assert(allowed_symmetries[idx].size() == characters[idx].size());
    for (int n_sym : allowed_symmetries[idx])
      assert(n_sym < n_symmetries_total_);

    n_symmetries_[name] = (int)allowed_symmetries[idx].size();
    allowed_symmetries_[name] = allowed_symmetries[idx];
    characters_[name] = characters[idx];

    // find out if characters are real
    is_real_[name] = true;
    for (auto c : characters[idx])
      if (std::abs(std::imag(c)) > 1e-12) {
        is_real_[name] = false;
        break;
      }
  }
}

std::vector<std::string> CharacterTable::names() const { return names_; }

SpaceGroup CharacterTable::little_group(std::string name) const {
  return space_group_.subgroup(allowed_symmetries_.find(name)->second);
}

int CharacterTable::n_symmetries(std::string name) const {
  return n_symmetries_.find(name)->second;
}

std::vector<int> CharacterTable::allowed_symmetries(std::string name) const {
  return allowed_symmetries_.find(name)->second;
}

complex CharacterTable::character(const std::string &name,
                                  const int &n_sym) const {
  return characters_.find(name)->second[n_sym];
}

std::vector<complex> CharacterTable::characters(const std::string &name) const {
  return characters_.find(name)->second;
}

// void Print(const CharacterTable& table)
// {
//   for (auto name : table.names())
// 	{
// 	  printf("[Representation]=%s\n", name.c_str());
// 	  if ( table.is_real(name) ) printf("REAL\n");
// 	  else printf("COMPLEX\n");
// 	  printf("[AllowedOps]=%d\n", table.n_symmetries(name));
// 	  for (int n_sym : table.allowed_symmetries(name))
// 	    printf("%d ", n_sym);
// 	  printf("\n");
// 	  for (int n_sym=0; n_sym < table.n_symmetries(name); ++n_sym)
// 	    printf("%f %f\n", std::real(table.character(name, n_sym)),
// 		   std::imag(table.character(name, n_sym)));
// 	}
// }

bool CharacterTable::is_real(const std::string &name) const {
  return is_real_.find(name)->second;
}

CharacterTable read_charactertable(std::string filename) {
  std::ifstream File(filename.c_str());
  if (File.fail()) {
    std::cerr << "Error in read_charactertable: Could not open "
              << "file with filename [" << filename << "] given. Abort."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  SpaceGroup space_group = read_spacegroup(filename);
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
  for (int i = 0; i < nreps; ++i) {
    File >> tobeparsed;
    while (tobeparsed.find("[Representation]") == std::string::npos)
      File >> tobeparsed;
    pos = tobeparsed.find('=');

    // Get name of representation
    std::string rep_name = tobeparsed.substr(pos + 1, std::string::npos);
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
        std::cerr << "Read Error in read_charactertable (I)" << std::endl;
        exit(EXIT_FAILURE);
      }
      allowed_ops.push_back(w);
    }

    if (!File.good()) {
      std::cerr << "Read Error in read_charactertable (II)" << std::endl;
      exit(EXIT_FAILURE);
    }

    // parse bloch factors
    for (int so = 0; so < n_allowed_ops; ++so) {
      double re, im;
      File >> re >> im;
      if (!File.good()) {
        std::cerr << "Read Error in read_charactertable (III)" << std::endl;
        exit(EXIT_FAILURE);
      }
      characters.push_back(re + complex(0, 1) * im);
    }
    names.push_back(rep_name);
    allowed_symmetries_arr.push_back(allowed_ops);
    characters_arr.push_back(characters);
  }

  return CharacterTable(space_group, names, allowed_symmetries_arr,
                        characters_arr);
}

} // namespace hydra
