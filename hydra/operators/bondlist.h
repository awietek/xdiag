#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "bond.h"

namespace hydra {

class BondList {
  using iterator_t = typename std::vector<Bond>::iterator;
  using const_iterator_t = typename std::vector<Bond>::const_iterator;

public:
  BondList() = default;
  explicit BondList(std::vector<Bond> const &bonds);

  void operator<<(Bond const &bond);

  bool coupling_defined(std::string name) const;
  void set_coupling(std::string name, complex cpl);
  complex get_coupling(std::string name) const;

  bool matrix_defined(std::string name) const;
  void set_matrix(std::string name, arma::cx_mat mat);
  void set_matrix(std::string name, arma::mat mat);
  arma::cx_mat get_matrix(std::string name) const;

  int size() const;
  void clear();
  int n_sites() const;

  BondList bonds_of_type(std::string type) const;
  BondList bonds_with_matrix() const;

  bool is_complex(double precision = 1e-12) const;
  bool is_real(double precision = 1e-12) const;

  inline iterator_t begin() { return bonds_.begin(); }
  inline iterator_t end() { return bonds_.end(); }
  inline const_iterator_t begin() const { return bonds_.begin(); }
  inline const_iterator_t end() const { return bonds_.end(); }
  inline const_iterator_t cbegin() const { return bonds_.cbegin(); }
  inline const_iterator_t cend() const { return bonds_.cend(); }

  friend BondList operator+(BondList const &, BondList const &);

private:
  std::vector<Bond> bonds_;
  std::map<std::string, complex> couplings_;
  std::map<std::string, arma::cx_mat> matrices_;
};

BondList operator+(BondList const &, BondList const &);
BondList read_bondlist(std::string filename);

} // namespace hydra
