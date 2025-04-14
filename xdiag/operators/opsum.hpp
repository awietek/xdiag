// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <map>
#include <string>
#include <vector>
#include <xdiag/operators/coupling.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

class OpSum {
public:
  using iterator_t = std::vector<std::pair<Coupling, Op>>::const_iterator;

  XDIAG_API OpSum() = default;
  XDIAG_API explicit OpSum(Op const &op);
  XDIAG_API OpSum(std::string coupling, Op const &op);
  XDIAG_API OpSum(double coupling, Op const &op);
  XDIAG_API OpSum(complex coupling, Op const &op);
  OpSum(Coupling const &coupling, Op const &op);

  XDIAG_API OpSum &operator=(Op const &op);

  XDIAG_API OpSum &operator*=(double scalar);
  XDIAG_API OpSum &operator*=(complex scalar);
  OpSum &operator*=(Scalar const &scalar);

  XDIAG_API OpSum &operator/=(double scalar);
  XDIAG_API OpSum &operator/=(complex scalar);
  OpSum &operator/=(Scalar const &scalar);

  XDIAG_API OpSum &operator+=(OpSum const &ops);
  XDIAG_API OpSum &operator+=(Op const &op);

  XDIAG_API OpSum operator+(OpSum const &ops) const;
  XDIAG_API OpSum operator+(Op const &op) const;

  XDIAG_API OpSum &operator-=(OpSum const &ops);
  XDIAG_API OpSum &operator-=(Op const &op);

  XDIAG_API OpSum operator-(OpSum const &ops) const;
  XDIAG_API OpSum operator-(Op const &op) const;

  XDIAG_API Scalar &operator[](std::string name);
  XDIAG_API Scalar const &operator[](std::string name) const;

  XDIAG_API bool operator==(OpSum const &rhs) const;
  XDIAG_API bool operator!=(OpSum const &rhs) const;

  std::vector<std::pair<Coupling, Op>> const &terms() const;
  std::vector<std::string> constants() const;
  XDIAG_API OpSum plain() const;
  XDIAG_API int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

private:
  std::vector<std::pair<Coupling, Op>> terms_;
  std::map<std::string, Scalar> constants_;
};

XDIAG_API std::vector<std::string> constants(OpSum const &ops);

// Creation
XDIAG_API OpSum operator*(double coupling, Op const &op);
XDIAG_API OpSum operator*(complex coupling, Op const &op);
XDIAG_API OpSum operator*(std::string coupling, Op const &op);

XDIAG_API OpSum operator*(Op const &op, double coupling);
XDIAG_API OpSum operator*(Op const &op, complex coupling);
XDIAG_API OpSum operator*(Op const &op, std::string coupling);

OpSum operator*(Scalar const &coupling, Op const &op);
OpSum operator*(Coupling const &coupling, Op const &op);
OpSum operator*(Op const &op, Scalar const &coupling);
OpSum operator*(Op const &op, Coupling const &coupling);

// scalar Multiplication
XDIAG_API OpSum operator*(double scalar, OpSum const &op);
XDIAG_API OpSum operator*(complex scalar, OpSum const &op);
OpSum operator*(Scalar const &coupling, OpSum const &op);

XDIAG_API OpSum operator*(OpSum const &op, double scalar);
XDIAG_API OpSum operator*(OpSum const &op, complex scalar);
OpSum operator*(OpSum const &op, Scalar const &coupling);

XDIAG_API OpSum operator/(OpSum const &op, double scalar);
XDIAG_API OpSum operator/(OpSum const &op, complex scalar);
OpSum operator/(OpSum const &op, Scalar const &coupling);

XDIAG_API std::ostream &operator<<(std::ostream &out, OpSum const &ops);
XDIAG_API std::string to_string(OpSum const &ops);

} // namespace xdiag
