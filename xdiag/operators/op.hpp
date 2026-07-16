// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>
#include <string>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/matrix.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Op represents a single elementary quantum operator, identified by:
//   - a type string (e.g. "Sz", "S+", "Hop", "Matrix")
//   - an optional ordered list of lattice sites it acts on
//   - an optional dense matrix for custom "Matrix" operators
//
// Known types and their site requirements are defined in logic/types.hpp.
// Use Monomial to form products of Ops, and OpSum for linear combinations.
class XDIAG_API Op {
public:
  Op() = default;
  explicit Op(std::string type); // site-free (e.g. "HubbardU")

  Op(std::string type, int64_t site); // single-site
  Op(std::string type,
     std::initializer_list<int64_t> const &sites); // multi-site
  Op(std::string type,
     std::vector<int64_t> const &sites); // multi-site
  Op(std::string type,
     arma::Col<int64_t> const &sites); // multi-site

  // Custom matrix operator on one or more sites (real or complex)
  Op(std::string type, int64_t site, arma::mat const &matrix);
  Op(std::string type, std::initializer_list<int64_t> const &sites,
     arma::mat const &matrix);
  Op(std::string type, std::vector<int64_t> const &sites,
     arma::mat const &matrix);
  Op(std::string type, arma::Col<int64_t> const &sites,
     arma::mat const &matrix);

  Op(std::string type, int64_t site, arma::cx_mat const &matrix);
  Op(std::string type, std::initializer_list<int64_t> const &sites,
     arma::cx_mat const &matrix);
  Op(std::string type, std::vector<int64_t> const &sites,
     arma::cx_mat const &matrix);
  Op(std::string type, arma::Col<int64_t> const &sites,
     arma::cx_mat const &matrix);
  Op(std::string type, int64_t site, Matrix const &matrix);
  Op(std::string type, std::vector<int64_t> const &sites, Matrix const &matrix);

  std::string type() const;              // operator type string
  int64_t size() const;                  // number of sites
  int64_t operator[](int64_t idx) const; // site at index idx
  std::vector<int64_t> const &sites() const;

  bool hassites() const;  // true if site list was provided
  bool hasmatrix() const; // true if a matrix was provided
  Matrix const &matrix() const;
  bool isreal() const;

  bool operator==(const Op &rhs) const;
  bool operator!=(const Op &rhs) const;
  bool operator<(const Op &rhs) const noexcept;

private:
  std::string type_;
  std::optional<std::vector<int64_t>> sites_;
  std::optional<Matrix> matrix_;
};

XDIAG_API bool isreal(Op const &op);
XDIAG_API std::ostream &operator<<(std::ostream &out, Op const &op);
XDIAG_API std::string to_string(Op const &op);

} // namespace xdiag
