#pragma once

#include <stdexcept>

namespace hydra {
class symmetry_error : public std::runtime_error {
public:
  symmetry_error(std::string const &what = "") : std::runtime_error(what) {}
};

} // namespace hydra
