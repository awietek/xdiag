#pragma once

#include <string>
#include <vector>

#include <hydra/common.h>
#include <lila/all.h>

namespace hydra {

class Representation {
public:
  Representation() = default;
  explicit Representation(std::vector<complex> const &characters);
  Representation(std::vector<complex> const &characters,
                 std::vector<int> const &allowed_symmetries);

  inline complex character(int idx) const { return characters_(idx); }
  std::vector<int> allowed_symmetries() const { return allowed_symmetries_; }

  Representation subgroup(std::vector<int> const &symmetry_numbers) const;

  inline idx_t size() const { return characters_.size(); }

  bool operator==(Representation const &rhs) const;
  bool operator!=(Representation const &rhs) const;

private:
  lila::Vector<complex> characters_;
  std::vector<int> allowed_symmetries_;
};

Representation read_represenation(std::string filename, std::string repname);

} // namespace hydra
