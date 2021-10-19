#pragma once

#include <map>
#include <string>
#include <vector>

#include <hydra/common.h>
#include <lila/all.h>

namespace hydra {

class Couplings {
  using iterator_t = typename std::map<std::string, complex>::iterator;
  using const_iterator_t =
      typename std::map<std::string, complex>::const_iterator;

public:
  Couplings() = default;

  explicit Couplings(const std::map<std::string, complex> &couplings)
      : couplings_(couplings) {}

  std::vector<std::string> couplings() const;

  complex operator[](const std::string &name) const {
    return couplings_.find(name)->second;
  }
  complex &operator[](const std::string &name) { return couplings_[name]; }

  bool defined(const std::string &name) const {
    return couplings_.find(name) != couplings_.end();
  }

  bool is_real(const std::string &name) const {
    return lila::close(std::abs(std::imag(couplings_.find(name)->second)), 0.);
  }

  bool all_real() const {
    for (auto name_val : couplings_) {
      if (!is_real(name_val.first))
        return false;
    }
    return true;
  }

  double real(const std::string &name) const {
    return std::real(couplings_.find(name)->second);
  }

  double imag(const std::string &name) const {
    return std::imag(couplings_.find(name)->second);
  }

  iterator_t begin() { return couplings_.begin(); }
  iterator_t end() { return couplings_.end(); }
  const_iterator_t begin() const { return couplings_.begin(); }
  const_iterator_t end() const { return couplings_.end(); }
  const_iterator_t cbegin() const { return couplings_.cbegin(); }
  const_iterator_t cend() const { return couplings_.cend(); }

  void clear() { couplings_.clear(); }

private:
  std::map<std::string, complex> couplings_;
};

Couplings read_couplings(std::string filename);

bool is_complex(Couplings const &cpls);
bool is_real(Couplings const &cpls);

} // namespace hydra
