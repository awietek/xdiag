#pragma once

#include <map>
#include <string>
#include <vector>

#include <hydra/common.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

class CharacterTable {
public:
  CharacterTable(std::vector<Representation> const &representations);

  std::vector<std::string> names() const { return names_; }
  int n_symmetries(std::string name) const { return n_symmetries_.at(name); }
  std::vector<int> symmetries(std::string name) const {
    return symmetries_.at(name);
  }
  std::vector<complex> characters(std::string name) const {
    return characters_.at(name);
  }
  bool is_real(std::string name) const { return is_real_.at(name); }
  std::vector<double> characters_real(std::string name) const {
    return characters_real_.at(name);
  }

private:
  std::vector<Representation> representations_;
  std::vector<std::string> names_;
  std::map<std::string, int> n_symmetries_;
  std::map<std::string, std::vector<int>> symmetries_;
  std::map<std::string, std::vector<complex>> characters_;
  std::map<std::string, bool> is_real_;
  std::map<std::string, std::vector<double>> characters_real_;
};

std::vector<Representation> read_represenations(std::string filename);
} // namespace hydra
