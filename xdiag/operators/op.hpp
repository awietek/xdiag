#pragma once

#include <vector>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/coupling.hpp>

namespace xdiag {

class Op {
public:
  Op() = default;
  Op(std::string type, Coupling const &coupling,
     std::vector<int64_t> const &sites);
  Op(std::string type, Coupling const &coupling, int64_t site);

  Op(std::string type, const char *coupling, std::vector<int64_t> const &sites);
  Op(std::string type, const char *coupling, int64_t site);
  Op(std::string type, std::string coupling, std::vector<int64_t> const &sites);
  Op(std::string type, std::string coupling, int64_t site);
  Op(std::string type, double coupling, std::vector<int64_t> const &sites);
  Op(std::string type, double coupling, int64_t site);
  Op(std::string type, complex coupling, std::vector<int64_t> const &sites);
  Op(std::string type, complex coupling, int64_t site);
  Op(std::string type, arma::mat const &coupling,
     std::vector<int64_t> const &sites);
  Op(std::string type, arma::mat const &coupling, int64_t site);
  Op(std::string type, arma::cx_mat const &coupling,
     std::vector<int64_t> const &sites);
  Op(std::string type, arma::cx_mat const &coupling, int64_t site);

  std::string type() const;
  Coupling const &coupling() const;

  int64_t size() const;
  int64_t operator[](int64_t idx) const;
  std::vector<int64_t> const &sites() const;

  bool isreal() const;
  bool ismatrix() const;
  bool isexplicit() const;

  bool operator==(const Op &rhs) const;
  bool operator!=(const Op &rhs) const;

private:
  std::string type_;
  Coupling coupling_;
  std::vector<int64_t> sites_;
};

bool sites_disjoint(Op const &op);
std::vector<int64_t> common_sites(Op const &b1, Op const &b2);

std::ostream &operator<<(std::ostream &out, Op const &op);
std::string to_string(Op const &op);

// helper for julia wrapper
class VectorOp {
public:
  VectorOp() = default;
  void push_back(Op const &op);
  std::vector<Op> vector() const;
private:
  std::vector<Op> v_;
};

} // namespace xdiag
