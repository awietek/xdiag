#include "opsum.hpp"

namespace xdiag {

OpSum::OpSum(Op const &op) : terms_({{Coupling(1.0), op}}) {}
OpSum::OpSum(Coupling const &cpl, Op const &op) : terms_({{cpl, op}}) {}

void OpSum::operator+=(OpSum const &ops) {
  terms_.insert(terms_.end(), ops.terms_.begin(), ops.terms_.end());
  constants_.insert(ops.constants_.begin(), ops.constants_.end());
}

OpSum OpSum::operator+(OpSum const &op) const {
  OpSum new_ops = *this;
  new_ops += op;
  return new_ops;
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
        ops += it->second * op;
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
