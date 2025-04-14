// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>
#include <xdiag/common.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

class Coupling {
public:
  using value_t = std::variant<std::string, Scalar>;

  XDIAG_API Coupling() = default;
  XDIAG_API explicit Coupling(std::string value);
  XDIAG_API explicit Coupling(const char *value);
  XDIAG_API explicit Coupling(double value);
  XDIAG_API explicit Coupling(complex value);
  XDIAG_API explicit Coupling(Scalar value);

  XDIAG_API bool operator==(Coupling const &rhs) const;
  XDIAG_API bool operator!=(Coupling const &rhs) const;

  XDIAG_API bool isscalar() const;
  XDIAG_API bool isstring() const;
  XDIAG_API Scalar scalar() const;
  XDIAG_API std::string string() const;

private:
  value_t value_;
};

XDIAG_API bool isscalar(Coupling const &c);
XDIAG_API bool isstring(Coupling const &c);
XDIAG_API Scalar scalar(Coupling const &c);
XDIAG_API std::string string(Coupling const &c);

XDIAG_API std::ostream &operator<<(std::ostream &out, Coupling const &v);
XDIAG_API std::string to_string(Coupling const &v);

} // namespace xdiag
