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
  XDIAG_API OpSum(Coupling const &cpl, Op const &op);
  XDIAG_API OpSum(std::string cpl, Op const &op);
  XDIAG_API OpSum(double cpl, Op const &op);
  XDIAG_API OpSum(complex cpl, Op const &op);

  XDIAG_API OpSum &operator=(Op const &op);

  XDIAG_API OpSum &operator*=(Scalar const &cpl);
  XDIAG_API OpSum &operator/=(Scalar const &cpl);
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
XDIAG_API OpSum operator*(Coupling const &cpl, Op const &op);
XDIAG_API OpSum operator*(Op const &op, Coupling const &cpl);
XDIAG_API OpSum operator*(Scalar const &cpl, OpSum const &op);
XDIAG_API OpSum operator*(OpSum const &op, Scalar const &cpl);
XDIAG_API OpSum operator/(OpSum const &op, Scalar const &cpl);

XDIAG_API std::ostream &operator<<(std::ostream &out, OpSum const &ops);
XDIAG_API std::string to_string(OpSum const &ops);

} // namespace xdiag
