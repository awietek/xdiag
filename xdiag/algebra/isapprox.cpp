// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isapprox.hpp"

#include <map>
#include <set>

#include <xdiag/math/matrix.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

bool isapprox(Op const &op1, Op const &op2, double rtol, double atol) {
  if (op1.type() != op2.type())
    return false;
  if (op1.hassites() != op2.hassites())
    return false;
  if (op1.hassites() && op1.sites() != op2.sites())
    return false;
  if (op1.hasmatrix() != op2.hasmatrix())
    return false;
  if (op1.hasmatrix())
    return xdiag::isapprox(op1.matrix(), op2.matrix(), rtol, atol);
  return true;
}

// Build a map {Monomial -> Scalar} summing coefficients for equal monomials,
// using ops.plain() to resolve named parameters.
static std::map<Monomial, Scalar> to_map(OpSum const &ops) {
  std::map<Monomial, Scalar> m;
  for (auto const &[c, mono] : ops.plain()) {
    auto it = m.find(mono);
    if (it == m.end())
      m[mono] = c.scalar();
    else
      it->second += c.scalar();
  }
  return m;
}

bool isapprox(OpSum const &ops1, OpSum const &ops2, double rtol,
              double atol) try {
  auto map1 = to_map(ops1);
  auto map2 = to_map(ops2);

  // Collect all monomials
  std::set<Monomial> all;
  for (auto const &[m, c] : map1)
    all.insert(m);
  for (auto const &[m, c] : map2)
    all.insert(m);

  for (auto const &m : all) {
    Scalar c1 = map1.count(m) ? map1.at(m) : Scalar(0.0);
    Scalar c2 = map2.count(m) ? map2.at(m) : Scalar(0.0);
    if (!xdiag::isapprox(c1, c2, rtol, atol))
      return false;
  }
  return true;
}
XDIAG_CATCH

std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        double rtol, double atol) try {
  auto map1 = to_map(ops1);
  auto map2 = to_map(ops2);

  std::set<Monomial> all;
  for (auto const &[m, c] : map1)
    all.insert(m);
  for (auto const &[m, c] : map2)
    all.insert(m);

  // Find lambda such that ops1 = lambda * ops2
  std::optional<Scalar> lambda;
  for (auto const &m : all) {
    Scalar c1 = map1.count(m) ? map1.at(m) : Scalar(0.0);
    Scalar c2 = map2.count(m) ? map2.at(m) : Scalar(0.0);

    bool z1 = xdiag::isapprox(c1, Scalar(0.0), rtol, atol);
    bool z2 = xdiag::isapprox(c2, Scalar(0.0), rtol, atol);

    if (z1 && z2)
      continue;
    if (z1 != z2)
      return std::nullopt;

    Scalar ratio = c1 / c2;
    if (!lambda) {
      lambda = ratio;
    } else {
      if (!xdiag::isapprox(*lambda, ratio, rtol, atol))
        return std::nullopt;
    }
  }
  return lambda;
}
XDIAG_CATCH

} // namespace xdiag
