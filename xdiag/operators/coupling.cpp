#include "coupling.hpp"

namespace xdiag {

Coupling(std::string value) : value_(value){};
Coupling(double value) : value_(Scalar(value)){};
Coupling(complex value) : value_(Scalar(value)){};
Coupling(Scalar value) : value_(value){};

bool isscalar() const { return std::holds_alternative<Scalar>(value_); }
bool isstring() const { return std::holds_alternative<std::string>(value_); }

Scalar scalar() const try {
  if (const Scalar *v = std::get_if<Scalar>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Coupling holding a value of type "
                "\"std::string\" to value of type \"Scalar\"");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::string string() const try {
  if (const std::string *v = std::get_if<std::string>(&value_)) {
    return *v;
  } else {
    XDIAG_THROW("Cannot convert Coupling holding a value of type "
                "\"Scalar\" to value of type \"std::string\"");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool operator==(Coupling const &rhs) const { return value_ == rhs.value_; }
bool operator!=(Coupling const &rhs) const { return !operator==(rhs); }

bool isscalar(Coupling const &c) { return c.isscalar(); }
bool isstring(Coupling const &c) { return c.isstring(); }
Scalar scalar(Coupling const &c) { return c.scalar(); }
std::string string(Coupling const &c) { return c.string(); }

std::ostream &operator<<(std::ostream &out, Coupling const &cpl) {
  std::visit([&](auto &&c) { out << c }, cpl);
}
std::string to_string(Coupling const &v) { return to_string_generic(v); }

} // namespace xdiag
