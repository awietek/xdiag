// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "op.hpp"

#include <string>

#include <xdiag/operators/types.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>
#include <xdiag/utils/xdiag_offset.hpp>

#ifndef XDIAG_DISABLE_COLOR
#include <extern/fmt/color.hpp>
#endif

namespace xdiag {

Op::Op(std::string type) : type_(type) {}
Op::Op(std::string type, int64_t site)
    : type_(type), sites_(std::vector<int64_t>{site}) {}
Op::Op(std::string type, std::vector<int64_t> const &sites)
    : type_(type), sites_(sites) {}
Op::Op(std::string type, arma::Col<int64_t> const &sites)
    : type_(type), sites_(std::vector<int64_t>(sites.begin(), sites.end())) {}
Op::Op(std::string type, std::initializer_list<int64_t> const &sites)
    : type_(type), sites_(std::vector<int64_t>(sites.begin(), sites.end())) {}

Op::Op(std::string type, int64_t site, arma::mat const &mat)
    : type_(type), sites_(std::vector<int64_t>{site}), matrix_(mat) {}
Op::Op(std::string type, std::vector<int64_t> const &sites,
       arma::mat const &mat)
    : type_(type), sites_(sites), matrix_(mat) {}
Op::Op(std::string type, arma::Col<int64_t> const &sites, arma::mat const &mat)
    : type_(type), sites_(std::vector<int64_t>(sites.begin(), sites.end())),
      matrix_(mat) {}
Op::Op(std::string type, std::initializer_list<int64_t> const &sites,
       arma::mat const &mat)
    : type_(type), sites_(std::vector<int64_t>(sites.begin(), sites.end())),
      matrix_(mat) {}

Op::Op(std::string type, int64_t site, arma::cx_mat const &mat)
    : type_(type), sites_(std::vector<int64_t>{site}), matrix_(mat) {}
Op::Op(std::string type, std::vector<int64_t> const &sites,
       arma::cx_mat const &mat)
    : type_(type), sites_(sites), matrix_(mat) {}
Op::Op(std::string type, arma::Col<int64_t> const &sites,
       arma::cx_mat const &mat)
    : type_(type), sites_(std::vector<int64_t>(sites.begin(), sites.end())),
      matrix_(mat) {}
Op::Op(std::string type, std::initializer_list<int64_t> const &sites,
       arma::cx_mat const &mat)
    : type_(type), sites_(std::vector<int64_t>(sites.begin(), sites.end())),
      matrix_(mat) {}

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
}
XDIAG_CATCH

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
}
XDIAG_CATCH

std::vector<int64_t> const &Op::sites() const try {
  if (hassites()) {
    return sites_.value();
  } else {
    XDIAG_THROW(
        fmt::format("Op does not have sites defined (type is \"{}\").", type_));
  }
}
XDIAG_CATCH

bool Op::isreal() const {
  if (hasmatrix()) {
    return matrix_->isreal();
  }
  return is_real_type(type_);
}

bool isreal(Op const &op) { return op.isreal(); }

bool Op::operator==(Op const &rhs) const {
  return (type_ == rhs.type_) && (sites_ == rhs.sites_) &&
         (matrix_ == rhs.matrix_);
}
bool Op::operator!=(Op const &rhs) const { return !operator==(rhs); }

bool Op::operator<(Op const &rhs) const noexcept {
  // 1. No-site ops come before ops with sites
  if (hassites() != rhs.hassites())
    return !hassites();
  // 2. Both have sites: fewer sites first, then lexicographic site values
  if (hassites()) {
    if (sites_->size() != rhs.sites_->size())
      return sites_->size() < rhs.sites_->size();
    if (*sites_ != *rhs.sites_)
      return *sites_ < *rhs.sites_;
  }
  // 3. Same site(s): compare by type string
  if (type_ != rhs.type_)
    return type_ < rhs.type_;
  // 4. Same type and sites: no-matrix before has-matrix, then by matrix value
  if (hasmatrix() != rhs.hasmatrix())
    return hasmatrix() < rhs.hasmatrix();
  if (hasmatrix())
    return matrix() < rhs.matrix();
  return false;
}

// Styled helpers — no-ops when XDIAG_DISABLE_COLOR is defined
static std::string styled_type(std::string const &s) {
#ifndef XDIAG_DISABLE_COLOR
  return fmt::format(fmt::emphasis::bold | fg(fmt::rgb(0xE06C75)), "{}", s);
#else
  return s;
#endif
}

static std::string styled_sites(std::string const &s) {
#ifndef XDIAG_DISABLE_COLOR
  return fmt::format(fg(fmt::rgb(0x61AFEF)), "{}", s);
#else
  return s;
#endif
}

std::ostream &operator<<(std::ostream &out, Op const &op) {
  out << styled_type(op.type());
  if (op.hassites()) {
    auto const &sites = op.sites();
    std::string s = "{";
    for (size_t i = 0; i < sites.size(); ++i) {
      if (i > 0)
        s += ",";
      s += std::to_string(sites[i] + XDIAG_OFFSET);
    }
    s += "}";
    out << styled_sites(s);
  }
  return out;
}

std::string to_string(Op const &op) { return utils::to_string_generic(op); }

} // namespace xdiag
