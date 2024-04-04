#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/io/toml/file_toml_handler.hpp>

namespace xdiag {

const std::vector<std::string> complex_bond_types = {"SCALARCHIRALITY"};

class Bond {
public:
  // Constructors with type name
  Bond() = default;
  Bond(std::string type, int64_t site);
  Bond(std::string type, std::vector<int64_t> const &sites);
  Bond(std::string type, double coupling, int64_t site);
  Bond(std::string type, double coupling, std::vector<int64_t> const &sites);
  Bond(std::string type, complex coupling, int64_t site);
  Bond(std::string type, complex coupling, std::vector<int64_t> const &sites);
  Bond(std::string type, std::string coupling_name, int64_t site);
  Bond(std::string type, std::string coupling_name,
       std::vector<int64_t> const &sites);

  // Constructors with bond matrix
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, int64_t site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, std::vector<int64_t> const &sites);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, double coupling, int64_t site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, double coupling,
       std::vector<int64_t> const &sites);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, complex coupling, int64_t site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, complex coupling,
       std::vector<int64_t> const &sites);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, std::string coupling_name,
       int64_t site);
  template <typename coeff_t>
  Bond(arma::Mat<coeff_t> const &matrix, std::string coupling_name,
       std::vector<int64_t> const &sites);

  explicit Bond(io::FileTomlHandler &&hdl);

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
  std::vector<int64_t> sites() const;
  int64_t site(int64_t j) const;
  int64_t size() const;
  int64_t operator[](int64_t j) const;

  bool iscomplex(double precision = 1e-12) const;
  bool isreal(double precision = 1e-12) const;

  bool operator==(const Bond &rhs) const;

private:
  std::string type_;
  arma::Mat<complex> matrix_;
  complex coupling_;
  std::string coupling_name_;
  std::vector<int64_t> sites_;
};

std::vector<int64_t> common_sites(Bond b1, Bond b2);
std::ostream &operator<<(std::ostream &out, const Bond &bond);

} // namespace xdiag
