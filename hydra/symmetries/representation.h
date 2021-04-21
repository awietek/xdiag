#pragma once

#include <string>
#include <vector>

#include <hydra/common.h>

namespace hydra {

struct Representation {
  std::string name;
  std::vector<complex> characters;
};

std::vector<Representation> read_represenations(std::string filename);
} // namespace hydra
