#pragma once

#include "extern/armadillo/armadillo"

#include <ostream>
#include <string>
#include <vector>

#include <hydra/common.h>

namespace hydra {

const std::vector<std::string> complex_bond_types = {"SCALARCHIRALITY"};

class Bond {
public:
  // Constructors with type name
  Bond() = default;
  Bond(std::string type, int site);
  Bond(std::string type, std::vector<int> const &sites);
  Bond(std::string type, complex coupling, int site);
  Bond(std::string type, complex coupling, std::vector<int> const &sites);
  Bond(std::string type, std::string coupling_name, int site);
  Bond(std::string type, std::string coupling_name,
       std::vector<int> const &sites);

  // Constructors with bond matrix
  template <typename coeff_t> Bond(arma::Mat<coeff_t> const &matrix, int site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, std::vector<int> const &sites);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, complex coupling, int site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, complex coupling,
       std::vector<int> const &sites);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, std::string coupling_name, int site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, std::string coupling_name,
       std::vector<int> const &sites);

  bool type_defined() const;
  bool matrix_defined() const;
  bool coupling_defined() const;
  bool coupling_named() const;
  bool sites_disjoint() const;

  std::string type() const;
  arma::cx_mat matrix() const;
  arma::mat matrix_real() const;

  template <typename coeff_t = complex>
  coeff_t coupling(double precision = 1e-12) const;

  std::string coupling_name() const;
  std::vector<int> sites() const;
  int site(int j) const;
  int size() const;
  int operator[](int j) const;

  bool is_complex(double precision = 1e-12) const;
  bool is_real(double precision = 1e-12) const;

  bool operator==(const Bond &rhs) const;

private:
  std::string type_;
  arma::Mat<complex> matrix_;
  complex coupling_;
  std::string coupling_name_;
  std::vector<int> sites_;
};

std::vector<int> common_sites(Bond b1, Bond b2);
std::ostream &operator<<(std::ostream &out, const Bond &bond);

} // namespace hydra
