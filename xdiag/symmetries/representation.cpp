// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representation.hpp"

#include <iomanip>

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

Representation::Representation(PermutationGroup const &group) try
    : Representation(group, arma::vec(group.size(), arma::fill::ones)) {}
XDIAG_CATCH

template <typename T>
static void check_characters(PermutationGroup const &group,
                             arma::Col<T> const &characters,
                             double precision = 1e-12) try {
  int64_t n = group.size();
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      T cij = characters[group.multiply(i, j)];
      T cicj = characters[i] * characters[j];
      if (std::abs(cij - cicj) > precision) {
        XDIAG_THROW(
            "Characters of Representation are not fulfilling multiplicative "
            "law with respect to the PermutationGroup, i.e. there exists group "
            "elements \"g\" and \"h\" such that c(g)*c(h) != c(gh)")
      }
    }
  }
}
XDIAG_CATCH

bool Representation::PermutationIrrep::operator==(
    PermutationIrrep const &rhs) const {
  return (group == rhs.group) && (characters == rhs.characters);
}
bool Representation::PermutationIrrep::operator!=(
    PermutationIrrep const &rhs) const {
  return !operator==(rhs);
}

Representation::Representation(PermutationGroup const &group,
                               Vector const &characters) try
    : type_("SitePermutation"), irrep_(PermutationIrrep{group, characters}) {
  if (group.size() != characters.size()) {
    XDIAG_THROW(
        "Size of PermutationGroup is not equal to number of characters given.");
  }

  if (characters.isreal()) {
    check_characters(group, characters.as<arma::vec>());
  } else {
    arma::cx_vec chars = characters.as<arma::cx_vec>();
    check_characters(group, chars);

    // If imaginary part is small, make it real
    double ni = arma::norm(arma::imag(chars));
    if (ni < 1e-14) {
      std::get<PermutationIrrep>(irrep_).characters = Vector(arma::real(chars));
    }
  }
}
XDIAG_CATCH

Representation::Representation(std::string type, int64_t charge) try
    : type_(type), irrep_(charge) {}
XDIAG_CATCH

bool Representation::operator==(Representation const &rhs) const {
  return (type_ == rhs.type_) && (irrep_ == rhs.irrep_);
}

bool Representation::operator!=(Representation const &rhs) const {
  return !operator==(rhs);
}

bool Representation::is_permutation() const {
  return std::holds_alternative<PermutationIrrep>(irrep_);
}
bool Representation::is_charge() const {
  return std::holds_alternative<int64_t>(irrep_);
}

bool Representation::isreal() const try {
  if (is_permutation()) {
    return std::get<PermutationIrrep>(irrep_).characters.isreal();
  } else {
    return true;
  }
}
XDIAG_CATCH

std::string Representation::type() const { return type_; }
PermutationGroup Representation::group() const try {
  if (is_permutation()) {
    return std::get<PermutationIrrep>(irrep_).group;
  } else {
    XDIAG_THROW(fmt::format(
        "Undefined PermutationGroup for Representation of type: \"{}\"",
        type_));
  }
}
XDIAG_CATCH

int64_t Representation::charge() const try {
  if (is_charge()) {
    return std::get<int64_t>(irrep_);
  } else {
    XDIAG_THROW(fmt::format(
        "Undefined charge for Representation of type: \"{}\"", type_));
  }
}
XDIAG_CATCH

Vector Representation::characters() const try {
  if (is_permutation()) {
    return std::get<PermutationIrrep>(irrep_).characters;
  } else {
    XDIAG_THROW(fmt::format(
        "Undefined characters for Representation of type: \"{}\"", type_));
  }
}
XDIAG_CATCH

bool isreal(Representation const &irrep) { return irrep.isreal(); }

bool isapprox(Representation const &r1, Representation const &r2,
              double tol) try {
  if (r1.type() != r2.type()) {
    return false;
  }
  if (r1.is_permutation() && r2.is_permutation()) {
    return (r1.group() == r2.group()) &&
           isapprox(r1.characters(), r2.characters(), tol, tol);
  } else if (r1.is_charge() && r2.is_charge()) {
    return r1.charge() == r2.charge();
  } else {
    return false;
  }
}
XDIAG_CATCH

Representation multiply(Representation const &r1,
                        Representation const &r2) try {
  if (r1.type() != r2.type()) {
    XDIAG_THROW(fmt::format("Cannot multiply two Representations of different "
                            "type: \"{}\" and \"{}\".",
                            r1.type(), r2.type()));
  }

  if (r1.is_permutation() && r2.is_permutation()) {
    // tensor product of site-permutation irreps: characters multiply
    if (r1.group() != r2.group()) {
      XDIAG_THROW("The two given Representations do not belong to the same "
                  "PermutationGroup.");
    }
    Vector c1 = r1.characters();
    Vector c2 = r2.characters();
    if (c1.isreal() && c2.isreal()) {
      arma::vec c = c1.as<arma::vec>() % c2.as<arma::vec>();
      return Representation(r1.group(), Vector(c));
    } else {
      // The Representation constructor narrows back to real if possible
      arma::cx_vec c = c1.as<arma::cx_vec>() % c2.as<arma::cx_vec>();
      return Representation(r1.group(), Vector(c));
    }
  } else if (r1.is_charge() && r2.is_charge()) {
    // tensor product of U(1)-type irreps: e^{i n theta} * e^{i m theta} adds
    // charges
    return Representation(r1.type(), r1.charge() + r2.charge());
  } else {
    XDIAG_THROW("Cannot multiply Representations of different kinds.");
  }
}
XDIAG_CATCH

Representation operator*(Representation const &r1,
                         Representation const &r2) try {
  return multiply(r1, r2);
}
XDIAG_CATCH

std::ostream &operator<<(std::ostream &out, Representation const &irrep) {
  out << "type      : " << irrep.type() << "\n";
  if (irrep.is_charge()) {
    out << "charge    : " << irrep.charge() << "\n";
    return out;
  }

  Vector characters = irrep.characters();
  out << "size      : " << characters.size() << "\n";
  out << "characters:\n";

  if (characters.is<arma::vec>()) {
    for (double c : characters.as<arma::vec>()) {
      out << std::setprecision(8) << c << "\n";
    }
  } else {
    for (complex c : characters.as<arma::cx_vec>()) {
      out << std::setprecision(8) << c << "\n";
    }
  }
  return out;
}
std::string to_string(Representation const &irrep) {
  return to_string_generic(irrep);
}

} // namespace xdiag
