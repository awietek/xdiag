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

  inline complex character(int idx) const { return characters_(idx); }
  inline idx_t size() const { return characters_.size(); }

  bool operator==(Representation const &rhs);
  bool operator!=(Representation const &rhs);

private:
  lila::Vector<complex> characters_;
};

Representation read_represenation(std::string filename, std::string repname);
  
} // namespace hydra
