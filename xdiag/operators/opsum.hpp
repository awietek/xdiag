#pragma once

#include <map>
#include <string>
#include <vector>
#include <xdiag/operators/coupling.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/scalar.hpp>

namespace xdiag {

class OpSum {
public:
  using iterator_t = std::vector<std::pair<Coupling, Op>>::const_iterator;

  XDIAG_API OpSum() = default;
  XDIAG_API OpSum(Op const &op);
  XDIAG_API OpSum(Coupling const &cpl, Op const &op);

  XDIAG_API void operator+=(OpSum const &ops);
  XDIAG_API OpSum operator+(OpSum const &ops) const;

  XDIAG_API Scalar &operator[](std::string name);
  XDIAG_API Scalar const &operator[](std::string name) const;

  bool operator==(OpSum const &rhs) const;
  bool operator!=(OpSum const &rhs) const;

  std::vector<std::string> constants() const;
  OpSum plain() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

private:
  std::vector<std::pair<Coupling, Op>> terms_;
  std::map<std::string, Scalar> constants_;
};

XDIAG_API std::vector<std::string> constants(OpSum const &ops);

XDIAG_API OpSum operator*(Coupling const &cpl, Op const &op);
XDIAG_API OpSum operator*(Op const &op, Coupling const &cpl);

XDIAG_API std::ostream &operator<<(std::ostream &out, OpSum const &ops);
XDIAG_API std::string to_string(OpSum const &ops);

} // namespace xdiag
