// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "op.hpp"

namespace xdiag {

Op::Op(std::string type) : type_(type) {}
Op::Op(std::string type, int64_t site)
    : type_(type), sites_(std::vector<int64_t>{site}) {}
Op::Op(std::string type, std::vector<int64_t> const &sites)
    : type_(type), sites_(sites) {}

Op::Op(std::string type, int64_t site, arma::mat const &mat)
    : type_(type), sites_(std::vector<int64_t>{site}), matrix_(mat) {}
Op::Op(std::string type, std::vector<int64_t> const &sites,
       arma::mat const &mat)
    : type_(type), sites_(sites), matrix_(mat) {}

Op::Op(std::string type, int64_t site, arma::cx_mat const &mat)
    : type_(type), sites_(std::vector<int64_t>{site}), matrix_(mat) {}
Op::Op(std::string type, std::vector<int64_t> const &sites,
       arma::cx_mat const &mat)
    : type_(type), sites_(sites), matrix_(mat) {}

Op::Op(std::string type, int64_t site, Matrix const &mat)
    : type_(type), sites_(std::vector<int64_t>{site}), matrix_(mat) {}
Op::Op(std::string type, std::vector<int64_t> const &sites, Matrix const &mat)
    : type_(type), sites_(sites), matrix_(mat) {}

std::string Op::type() const { return type_; }
bool Op::hasmatrix() const { return matrix_.has_value(); }
bool Op::hassites() const { return sites_.has_value(); }

Matrix const &Op::matrix() const try {
  if (hasmatrix()) {
    return matrix_.value();
  } else {
    XDIAG_THROW(fmt::format(
        "Op does not have a matrix defined (type is \"{}\").", type_));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
int64_t Op::size() const { return sites_ ? sites_->size() : 0; }
int64_t Op::operator[](int64_t idx) const try {
  if (sites_) {
    try {
      return sites_->at(idx);
    } catch (std::out_of_range const &exc) {
      XDIAG_THROW(fmt::format(
          "Site index \"{}\" out of range for Op of size \"{}\"", idx, size()));
    }
  } else {
    XDIAG_THROW("Cannot access site of Op, since Op has no sites defined.");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::vector<int64_t> const &Op::sites() const try {
  if (hassites()) {
    return sites_.value();
  } else {
    XDIAG_THROW(
        fmt::format("Op does not have sites defined (type is \"{}\").", type_));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool Op::operator==(Op const &rhs) const {
  return (type_ == rhs.type_) && (sites_ == rhs.sites_) &&
         (matrix_ == rhs.matrix_);
}
bool Op::operator!=(Op const &rhs) const { return !operator==(rhs); }

std::ostream &operator<<(std::ostream &out, Op const &op) {
  out << op.type();
  if (op.hassites()) {
    out << " ";
    for (int64_t site : op.sites()) {
      out << site + XDIAG_OFFSET << " ";
    }
  }
  if (op.hasmatrix()) {
    out << "\n" << op.matrix();
  }
  return out;
}

std::string to_string(Op const &op) { return to_string_generic(op); }

void VectorOp::push_back(Op const &op) { v_.push_back(op); }
std::vector<Op> VectorOp::vector() const { return v_; }

} // namespace xdiag
