// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representation_set.hpp"

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

RepresentationSet::RepresentationSet(
    std::vector<Representation> const &representations) try
    : representations_(representations) {
  // Behave like a set keyed by type: reject duplicate types.
  int64_t n = (int64_t)representations_.size();
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = i + 1; j < n; ++j) {
      if (representations_[i].type() == representations_[j].type()) {
        XDIAG_THROW(fmt::format("RepresentationSet contains more than one "
                                "Representation of type \"{}\".",
                                representations_[i].type()));
      }
    }
  }
}
XDIAG_CATCH

RepresentationSet::RepresentationSet(
    std::initializer_list<Representation> representations) try
    : RepresentationSet(std::vector<Representation>(representations)) {}
XDIAG_CATCH

bool RepresentationSet::has_type(std::string type) const {
  for (Representation const &rep : representations_) {
    if (rep.type() == type) {
      return true;
    }
  }
  return false;
}

std::optional<PermutationGroup>
RepresentationSet::group(std::string type) const {
  for (Representation const &rep : representations_) {
    if ((rep.type() == type) && rep.is_permutation()) {
      return rep.group();
    }
  }
  return std::nullopt;
}

std::optional<int64_t> RepresentationSet::charge(std::string type) const {
  for (Representation const &rep : representations_) {
    if ((rep.type() == type) && rep.is_charge()) {
      return rep.charge();
    }
  }
  return std::nullopt;
}

std::optional<Vector> RepresentationSet::characters(std::string type) const {
  for (Representation const &rep : representations_) {
    if ((rep.type() == type) && rep.is_permutation()) {
      return rep.characters();
    }
  }
  return std::nullopt;
}

bool RepresentationSet::operator==(RepresentationSet const &rhs) const {
  if (representations_.size() != rhs.representations_.size()) {
    return false;
  }
  // order-independent: every Representation must have an equal counterpart of
  // the same type in rhs
  for (Representation const &rep : representations_) {
    bool found = false;
    for (Representation const &rep_rhs : rhs.representations_) {
      if (rep.type() == rep_rhs.type()) {
        if (rep != rep_rhs) {
          return false;
        }
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

bool RepresentationSet::operator!=(RepresentationSet const &rhs) const {
  return !operator==(rhs);
}

bool RepresentationSet::isreal() const try {
  for (Representation const &rep : representations_) {
    if (!rep.isreal()) {
      return false;
    }
  }
  return true;
}
XDIAG_CATCH

int64_t RepresentationSet::size() const {
  return (int64_t)representations_.size();
}

RepresentationSet::const_iterator RepresentationSet::begin() const {
  return representations_.begin();
}
RepresentationSet::const_iterator RepresentationSet::end() const {
  return representations_.end();
}

bool RepresentationSet::isapprox(RepresentationSet const &rhs, double tol) const
    try {
  if (representations_.size() != rhs.representations_.size()) {
    return false;
  }
  for (Representation const &rep : representations_) {
    bool found = false;
    for (Representation const &rep_rhs : rhs.representations_) {
      if (rep.type() == rep_rhs.type()) {
        if (!xdiag::isapprox(rep, rep_rhs, tol)) {
          return false;
        }
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}
XDIAG_CATCH

RepresentationSet
RepresentationSet::multiply(RepresentationSet const &rhs) const try {
  if (representations_.size() != rhs.representations_.size()) {
    XDIAG_THROW("Cannot multiply two RepresentationSets holding a different "
                "number of Representations.");
  }
  std::vector<Representation> result;
  result.reserve(representations_.size());
  for (Representation const &rep : representations_) {
    bool found = false;
    for (Representation const &rep_rhs : rhs.representations_) {
      if (rep.type() == rep_rhs.type()) {
        result.push_back(xdiag::multiply(rep, rep_rhs));
        found = true;
        break;
      }
    }
    if (!found) {
      XDIAG_THROW(
          fmt::format("Cannot multiply RepresentationSets: no "
                      "Representation of type \"{}\" in the second set.",
                      rep.type()));
    }
  }
  return RepresentationSet(result);
}
XDIAG_CATCH

bool isreal(RepresentationSet const &set) { return set.isreal(); }

bool isapprox(RepresentationSet const &s1, RepresentationSet const &s2,
              double tol) try {
  return s1.isapprox(s2, tol);
}
XDIAG_CATCH

RepresentationSet multiply(RepresentationSet const &s1,
                           RepresentationSet const &s2) try {
  return s1.multiply(s2);
}
XDIAG_CATCH

RepresentationSet operator*(RepresentationSet const &s1,
                            RepresentationSet const &s2) try {
  return multiply(s1, s2);
}
XDIAG_CATCH

} // namespace xdiag
