#pragma once

#include <string>
#include <variant>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

using coupling_t =
    std::variant<std::string, double, complex, arma::mat, arma::cx_mat>;

class Bond {
public:
  // Constructors with type name
  Bond() = default;
  Bond(std::string type, coupling_t coupling,
       std::vector<int64_t> const &sites);
  Bond(std::string type, coupling_t coupling, int64_t site);

  std::string type() const;
  template <typename T> bool coupling_is() const;
  template <typename T> T coupling() const;
  coupling_t coupling() const;
  std::string coupling_type_string() const;

  int64_t size() const;
  int64_t operator[](int64_t idx) const;
  std::vector<int64_t> const &sites() const;

  bool isexplicit() const;
  bool isreal() const;
  bool ismatrix() const;

  bool operator==(const Bond &rhs) const;
  bool operator!=(const Bond &rhs) const;

private:
  std::string type_;
  coupling_t coupling_;
  std::vector<int64_t> sites_;
};

bool sites_disjoint(Bond const &bond);
std::vector<int64_t> common_sites(Bond const &b1, Bond const &b2);
std::ostream &operator<<(std::ostream &out, Bond const &bond);


} // namespace xdiag
