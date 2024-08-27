#include "op.hpp"

#include <algorithm>
#include <set>

#include <xdiag/common.hpp>
#include <xdiag/utils/close.hpp>

namespace xdiag {

Op::Op(std::string type, Coupling const &coupling,
       std::vector<int64_t> const &sites) try
    : type_(type), coupling_(coupling), sites_(sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Op::Op(std::string type, Coupling const &coupling, int64_t site) try
    : type_(type), coupling_(coupling), sites_({site}) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Op::Op(std::string type, const char *coupling,
       std::vector<int64_t> const &sites) try
    : Op(type, Coupling(std::string(coupling)), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, const char *coupling, int64_t site) try
    : Op(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Op::Op(std::string type, std::string coupling,
       std::vector<int64_t> const &sites) try
    : Op(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, std::string coupling, int64_t site) try
    : Op(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, double coupling, std::vector<int64_t> const &sites) try
    : Op(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, double coupling, int64_t site) try
    : Op(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, complex coupling,
       std::vector<int64_t> const &sites) try
    : Op(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, complex coupling, int64_t site) try
    : Op(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Op::Op(std::string type, arma::mat const &coupling,
       std::vector<int64_t> const &sites) try
    : Op(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, arma::mat const &coupling, int64_t site) try
    : Op(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Op::Op(std::string type, arma::cx_mat const &coupling,
       std::vector<int64_t> const &sites) try
    : Op(type, Coupling(coupling), sites) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
Op::Op(std::string type, arma::cx_mat const &coupling, int64_t site) try
    : Op(type, coupling, std::vector<int64_t>({site})) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::string Op::type() const { return type_; }

Coupling const &Op::coupling() const { return coupling_; }

int64_t Op::size() const { return sites_.size(); }
int64_t Op::operator[](int64_t idx) const try {
  return sites_.at(idx);
} catch (std::out_of_range const &exc) {
  XDIAG_THROW("Site index out of range for Op");
}
std::vector<int64_t> const &Op::sites() const { return sites_; }

bool Op::isexplicit() const try {
  return coupling_.isexplicit();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Op::isreal() const try {
  if (coupling_.is<std::string>()) {
    XDIAG_THROW("Cannot determine whether coupling is real, since coupling is "
                "set to a std::string value");
  } else if (type_ == "SCALARCHIRALITY") {
    return false;
  } else {
    return coupling_.isreal();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Op::ismatrix() const { return coupling_.ismatrix(); }

bool Op::operator==(Op const &rhs) const try {
  return (type_ == rhs.type_) && (coupling_ == rhs.coupling_) &&
         (sites_ == rhs.sites_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Op::operator!=(Op const &rhs) const { return !operator==(rhs); }

bool sites_disjoint(Op const &op) {
  auto const &sites = op.sites();
  auto set = std::set<int64_t>(sites.begin(), sites.end());
  return set.size() == sites.size();
}

std::vector<int64_t> common_sites(Op const &b1, Op const &b2) {
  std::vector<int64_t> s1 = b1.sites();
  std::vector<int64_t> s2 = b2.sites();
  std::vector<int64_t> s12;
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(s12));
  return s12;
}

std::ostream &operator<<(std::ostream &out, Op const &op) {
  if (op.ismatrix()) {
    out << op.type() << " ";
    for (int64_t site : op.sites()) {
      out << site << " ";
    }
    out << "\n";
    out << op.coupling();
  } else {
    out << op.type() << " " << op.coupling() << " ";
    for (int64_t site : op.sites()) {
      out << site << " ";
    }
    out << "\n";
  }
  return out;
}

std::string to_string(Op const &op) { return to_string_generic(op); }

void VectorOp::push_back(Op const &op) { v_.push_back(op); }
std::vector<Op> VectorOp::vector() const { return v_; }

} // namespace xdiag
