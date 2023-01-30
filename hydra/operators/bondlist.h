#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist_handler.h>

namespace hydra {

class BondList {
  using iterator_t = typename std::vector<Bond>::iterator;
  using const_iterator_t = typename std::vector<Bond>::const_iterator;

public:
  BondList() = default;
  explicit BondList(std::vector<Bond> const &bonds);

  int size() const;
  void clear();
  int n_sites() const;

  bool coupling_defined(std::string name) const;
  void set_coupling(std::string name, complex cpl);
  template <typename coeff_t = complex>
  coeff_t coupling(std::string name, double precision = 1e-12) const;
  std::map<std::string, complex> const &couplings() const;

  bool matrix_defined(std::string name) const;
  void set_matrix(std::string name, arma::cx_mat mat);
  void set_matrix(std::string name, arma::mat mat);
  arma::cx_mat matrix(std::string name) const;
  std::map<std::string, arma::cx_mat> const &matrices() const;

  BondListHandler operator[](std::string name);
  bool operator==(BondList const &other) const;
  bool operator!=(BondList const &other) const;
  BondList operator+(BondList const & other);
  void operator<<(Bond const &bond);

  BondList bonds_of_type(std::string type) const;
  BondList bonds_with_matrix() const;

  bool is_complex(double precision = 1e-12) const;
  bool is_real(double precision = 1e-12) const;

  iterator_t begin();
  iterator_t end();
  const_iterator_t begin() const;
  const_iterator_t end() const;
  const_iterator_t cbegin() const;
  const_iterator_t cend() const;

private:
  std::vector<Bond> bonds_;
  std::map<std::string, complex> couplings_;
  std::map<std::string, arma::cx_mat> matrices_;
};

BondList read_bondlist(std::string filename);

} // namespace hydra
