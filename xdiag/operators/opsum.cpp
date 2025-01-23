#include "opsum.hpp"

#include <xdiag/operators/logic/order.hpp>

namespace xdiag {

OpSum::OpSum(Op const &op) : terms_({{Coupling(1.0), op}}) {}
OpSum::OpSum(Coupling const &cpl, Op const &op) : terms_({{cpl, op}}) {}

OpSum &OpSum::operator=(Op const &op) {
  terms_ = std::vector<std::pair<Coupling, Op>>{{Coupling(1.0), op}};
  constants_.clear();
  return *this;
}

OpSum &OpSum::operator*=(Scalar const &cpl) try {
  for (auto &[c, op] : terms_) {
    if (c.isscalar()) {
      c = c.scalar() * cpl;
    } else {
      auto it = constants_.find(c.string());
      if (it != constants_.end()) {
        c = it->second * cpl;
      } else {
        XDIAG_THROW(fmt::format("Cannot multiply OpSum with a Scalar. Coupling "
                                "given by string \"{}\" has not been defined",
                                c.string()));
      }
    }
  }
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum &OpSum::operator/=(Scalar const &cpl) try {
  return operator*=(Scalar(1.0) / cpl);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum &OpSum::operator+=(OpSum const &ops) try {
  terms_.insert(terms_.end(), ops.terms_.begin(), ops.terms_.end());

  // Check common constants have the same numerical value
  std::map<std::string, Scalar> new_constants;
  std::map<std::string, Scalar> const &cs1 = constants_;
  std::map<std::string, Scalar> const &cs2 = ops.constants_;

  // Add keys from cs1 and common keys
  for (auto &entry : cs1) {
    std::string key = entry.first;
    Scalar value1 = entry.second;

    // key not a member of cs2
    if (cs2.find(key) == cs2.end()) {
      new_constants[key] = value1;
    } else { // key is shared
      Scalar value2 = cs2.at(key);
      if (value1 != value2) {
        XDIAG_THROW("Conflicting values for coupling constant \"{}\"");
      } else {
        new_constants[key] = value1;
      }
    }
  }

  // Add keys from cs2 exclusively
  for (auto &entry : cs2) {
    std::string key = entry.first;
    Scalar value = entry.second;

    // key not a member of cs2
    if (cs1.find(key) == cs1.end()) {
      new_constants[key] = value;
    }
  }
  constants_ = new_constants;
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum &OpSum::operator+=(Op const &op) try {
  return operator+=(OpSum(op));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum OpSum::operator+(OpSum const &op) const try {
  OpSum new_ops = *this;
  new_ops += op;
  return new_ops;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
OpSum OpSum::operator+(Op const &op) const try {
  return operator+(OpSum(op));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum &OpSum::operator-=(OpSum const &ops) try {
  auto nops = ops * Scalar(-1.0);
  return operator+=(nops);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
OpSum &OpSum::operator-=(Op const &op) try {
  return operator-=(OpSum(op));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum OpSum::operator-(OpSum const &ops) const try {
  OpSum new_ops = *this;
  new_ops -= ops;
  return new_ops;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum OpSum::operator-(Op const &op) const try {
  return operator-(OpSum(op));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t OpSum::size() const { return terms_.size(); }

std::vector<std::pair<Coupling, Op>> const &OpSum::terms() const {
  return terms_;
}
std::vector<std::string> OpSum::constants() const {
  std::vector<std::string> names;
  for (auto const &c : constants_) {
    names.push_back(c.first);
  }
  return names;
}

Scalar &OpSum::operator[](std::string name) { return constants_[name]; }
Scalar const &OpSum::operator[](std::string name) const {
  return constants_.at(name);
}

OpSum OpSum::plain() const try {
  OpSum ops;
  for (auto const &[cpl, op] : terms_) {
    if (cpl.isscalar()) {
      ops += cpl * op;
    } else {
      auto it = constants_.find(cpl.string());
      if (it != constants_.end()) {
        ops += Coupling(it->second) * op;
      } else {
        XDIAG_THROW(fmt::format("Cannot make OpSum plain, i.e. replace string "
                                "couplings with scalars. Coupling given by "
                                "string \"{}\" has not been defined",
                                cpl.string()));
      }
    }
  }
  return ops;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool OpSum::operator==(OpSum const &rhs) const {
  return (terms_ == rhs.terms_) && (constants_ == rhs.constants_);
}
bool OpSum::operator!=(OpSum const &rhs) const { return !operator==(rhs); }

OpSum::iterator_t OpSum::begin() const { return terms_.begin(); }
OpSum::iterator_t OpSum::end() const { return terms_.end(); }

OpSum operator*(std::string cpl, Op const &op) {
  return OpSum(Coupling(cpl), op);
}
OpSum operator*(Coupling const &cpl, Op const &op) { return OpSum(cpl, op); }
OpSum operator*(Op const &op, Coupling const &cpl) { return cpl * op; }

OpSum operator*(Scalar const &cpl, OpSum const &op) {
  auto newop = op;
  newop *= cpl;
  return newop;
}
OpSum operator*(OpSum const &op, Scalar const &cpl) { return cpl * op; }
OpSum operator/(OpSum const &op, Scalar const &cpl) {
  return op * (Scalar(1.0) / cpl);
}

std::ostream &operator<<(std::ostream &out, OpSum const &ops) {
  out << "Interactions:\n";
  out << "-------------\n";

  for (auto [cpl, op] : ops) {
    out << cpl << " " << op << "\n";
  }
  if (ops.constants().size() > 0) {
    out << "\nConstants:\n";
    out << "----------\n";

    for (auto name : ops.constants()) {
      out << name << ": " << ops[name] << "\n";
    }
  }
  return out;
}

std::string to_string(OpSum const &ops) { return to_string_generic(ops); }

} // namespace xdiag
