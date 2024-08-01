#include "opsum.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

namespace xdiag {

OpSum::OpSum(std::vector<Op> const &ops) : ops_(ops) {}
int64_t OpSum::size() const { return ops_.size(); }
Coupling &OpSum::operator[](std::string name) { return couplings_[name]; }
Coupling const &OpSum::operator[](std::string name) const {
  return couplings_.at(name);
}

bool OpSum::defined(std::string name) const {
  return couplings_.find(name) != couplings_.end();
}

std::vector<std::string> OpSum::couplings() const {
  std::vector<std::string> names;
  for (auto const &cpl : couplings_) {
    names.push_back(cpl.first);
  }
  return names;
}

bool OpSum::isreal() const try {
  OpSum ops_explicit = make_explicit(*this);
  return std::all_of(ops_explicit.begin(), ops_explicit.end(),
                     [](Op const &b) { return b.isreal(); });
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool OpSum::isexplicit() const try {
  return std::all_of(ops_.begin(), ops_.end(),
                     [](Op const &b) { return b.isexplicit(); });
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void OpSum::push_back(Op const &op) { ops_.push_back(op); }
void OpSum::operator+=(Op const &op) { push_back(op); }
void OpSum::operator+=(OpSum const &ops) try {
  // Add ops
  for (auto const &op : ops.ops_) {
    ops_.push_back(op);
  }

  // Add possible couplings
  for (auto it = ops.couplings_.begin(); it != ops.couplings_.end(); it++) {
    std::string name = it->first;
    Coupling cpl = it->second;

    // if name does not exist yet, add coupling
    if (couplings_.find(name) == couplings_.end()) {
      couplings_[name] = cpl;
    }
    // If it exists, throw an error if the coupling disagrees
    else {
      if (couplings_[name] != cpl) {
        XDIAG_THROW(fmt::format("Conflicting values of coupling \"{}\" found "
                                "when trying to add two OpSums.",
                                name));
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum OpSum::operator+(Op const &op) const {
  OpSum new_ops = *this;
  new_ops += op;
  return new_ops;
}

OpSum OpSum::operator+(OpSum const &ops) const try {
  OpSum new_ops = *this;
  new_ops += ops;
  return new_ops;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool OpSum::operator==(OpSum const &other) const {
  if (ops_ != other.ops_) {
    return false;
  }
  auto c1s = couplings();
  auto c2s = other.couplings();
  if (c1s != c2s) {
    return false;
  }
  for (std::string c : c1s) {
    if (couplings_.at(c) != other.couplings_.at(c)) {
      return false;
    }
  }
  return true;
}
bool OpSum::operator!=(OpSum const &other) const { return !operator==(other); }

OpSum::iterator_t OpSum::begin() const { return ops_.begin(); }
OpSum::iterator_t OpSum::end() const { return ops_.end(); }

OpSum make_explicit(OpSum const &ops) try {
  OpSum explicit_ops;
  for (auto const &op : ops) {
    if (op.isexplicit()) {
      explicit_ops += op;
    } else {
      std::string name = op.coupling().as<std::string>();
      if (ops.defined(name)) {
        std::string type = op.type();
        Coupling cpl = ops[name];
        if (!cpl.isexplicit()) {
          XDIAG_THROW(fmt::format("Unable to make OpSum explicit: coupling "
                                  "\"{}\" is yet a string",
                                  name));
        }
        explicit_ops += Op(type, cpl, op.sites());
      } else {
        XDIAG_THROW(fmt::format(
            "Unable to make OpSum explicit: coupling \"{}\" is not defined",
            name));
      }
    }
  }
  return explicit_ops;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum ops_of_type(std::string type, OpSum const &ops) {
  OpSum new_ops;
  for (auto const &op : ops) {
    if (op.type() == type) {
      new_ops += op;
    }
  }
  return new_ops;
}

OpSum read_opsum(std::string filename) {
  std::vector<Op> ops;

  // Open file and handle error
  std::ifstream File(filename.c_str());
  if (File.fail()) {
    std::cerr << "Error in read_oplist: "
              << "Could not open file with filename [" << filename
              << "] given. Abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Advance to interaction lines
  std::string tobeparsed;
  getline(File, tobeparsed);
  while ((tobeparsed.find("[Interactions]") == std::string::npos) &&
         (tobeparsed.find("[interactions]") == std::string::npos))
    getline(File, tobeparsed);

  // read lines until '[' is found or else until EOF
  while (std::getline(File, tobeparsed)) {
    if ((tobeparsed.find('[') != std::string::npos))
      break;

    std::string type, coupling;
    std::vector<int64_t> sites;
    std::stringstream stream(tobeparsed);
    stream >> type;
    stream >> coupling;

    // Just parse sites
    int64_t n;
    while (stream >> n)
      sites.push_back(n);
    ops.push_back(Op(type, coupling, sites));

    if (!File.good())
      break;
  }

  return OpSum(ops);
}

std::ostream &operator<<(std::ostream &out, OpSum const &ops) {
  out << "Interactions:\n";
  out << "-------------\n";

  for (auto op : ops) {
    out << op;
  }
  if (ops.couplings().size() > 0) {
    out << "Couplings:\n";
    out << "----------\n";

    for (auto name : ops.couplings()) {
      Coupling cpl = ops[name];
      out << name << ": " << cpl << "\n";
    }
  }
  return out;
}

std::string to_string(OpSum const &ops) { return to_string_generic(ops); }

} // namespace xdiag
