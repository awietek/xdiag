// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coupling.hpp"

namespace xdiag {

Coupling::Coupling(std::string value) : value_(value) {}
Coupling::Coupling(const char *value) : value_(std::string(value)) {}
Coupling::Coupling(double value) : value_(Scalar(value)) {}
Coupling::Coupling(complex value) : value_(Scalar(value)) {}
Coupling::Coupling(Scalar value) : value_(value) {}

bool Coupling::isscalar() const {
  return std::holds_alternative<Scalar>(value_);
}
bool Coupling::isstring() const {
  return std::holds_alternative<std::string>(value_);
}

Scalar Coupling::scalar() const try {
  if (const Scalar *v = std::get_if<Scalar>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Coupling holding a value of type "
                "\"std::string\" to value of type \"Scalar\"");
  }
}
XDIAG_CATCH

std::string Coupling::string() const try {
  if (const std::string *v = std::get_if<std::string>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Coupling holding a value of type "
                "\"Scalar\" to value of type \"std::string\"");
  }
}
XDIAG_CATCH

bool Coupling::operator==(Coupling const &rhs) const {
  return value_ == rhs.value_;
}
bool Coupling::operator!=(Coupling const &rhs) const {
  return !operator==(rhs);
}

bool isscalar(Coupling const &c) { return c.isscalar(); }
bool isstring(Coupling const &c) { return c.isstring(); }
Scalar scalar(Coupling const &c) { return c.scalar(); }
std::string string(Coupling const &c) { return c.string(); }

std::ostream &operator<<(std::ostream &out, Coupling const &cpl) {
  if (cpl.isscalar()) {
    out << cpl.scalar();
  } else {
    out << cpl.string();
  }
  return out;
}
std::string to_string(Coupling const &v) { return to_string_generic(v); }

} // namespace xdiag
