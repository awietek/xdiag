#pragma once

#include <vector>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/coupling.hpp>

namespace xdiag {

class Bond {
public:
  // Constructors with type name
  Bond() = default;
  Bond(std::string type, Coupling coupling, std::vector<int64_t> const &sites);
  Bond(std::string type, Coupling coupling, int64_t site);

  Bond(std::string type, const char *coupling,
       std::vector<int64_t> const &sites);
  Bond(std::string type, const char *coupling, int64_t site);
  Bond(std::string type, std::string coupling,
       std::vector<int64_t> const &sites);
  Bond(std::string type, std::string coupling, int64_t site);
  Bond(std::string type, double coupling, std::vector<int64_t> const &sites);
  Bond(std::string type, double coupling, int64_t site);
  Bond(std::string type, complex coupling, std::vector<int64_t> const &sites);
  Bond(std::string type, complex coupling, int64_t site);
  Bond(std::string type, arma::mat const &coupling,
       std::vector<int64_t> const &sites);
  Bond(std::string type, arma::mat const &coupling, int64_t site);
  Bond(std::string type, arma::cx_mat const &coupling,
       std::vector<int64_t> const &sites);
  Bond(std::string type, arma::cx_mat const &coupling, int64_t site);

  std::string type() const;
  Coupling const &coupling() const;

  int64_t size() const;
  int64_t operator[](int64_t idx) const;
  std::vector<int64_t> const &sites() const;

  bool isreal() const;
  bool ismatrix() const;
  bool isexplicit() const;

  bool operator==(const Bond &rhs) const;
  bool operator!=(const Bond &rhs) const;

private:
  std::string type_;
  Coupling coupling_;
  std::vector<int64_t> sites_;
};

bool sites_disjoint(Bond const &bond);
std::vector<int64_t> common_sites(Bond const &b1, Bond const &b2);
std::ostream &operator<<(std::ostream &out, Bond const &bond);

} // namespace xdiag
