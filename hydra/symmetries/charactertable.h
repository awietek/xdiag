#pragma once

#include <map>
#include <string>
#include <vector>

#include <hydra/common.h>
#include <hydra/symmetries/spacegroup.h>

namespace hydra {
class CharacterTable {

public:
  CharacterTable(SpaceGroup const& space_group,
                 std::vector<std::string> const& names,
                 std::vector<std::vector<int>> const& allowed_symmetries,
                 std::vector<std::vector<complex>> const& characters);

  SpaceGroup space_group() const { return space_group_; }
  std::vector<std::string> names() const;
  SpaceGroup little_group(std::string name) const;
  int n_symmetries(std::string name) const;
  std::vector<int> allowed_symmetries(std::string name) const;
  complex character(std::string name, int n_sym) const;
  std::vector<complex> characters(std::string name) const;
  std::vector<double> characters_real(std::string name) const;
  bool is_real(std::string name) const;

private:
  SpaceGroup space_group_;
  std::vector<std::string> names_;
  int n_symmetries_total_;
  std::map<std::string, int> n_symmetries_;
  std::map<std::string, std::vector<int>> allowed_symmetries_;
  std::map<std::string, std::vector<complex>> characters_;
  std::map<std::string, bool> is_real_;
};

CharacterTable read_charactertable(std::string filename);
} // namespace hydra
