// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <variant>

#include <xdiag/math/vector.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class Representation {
public:
  // A 1-D irrep of a finite group of site permutations: one character per
  // group element.
  struct PermutationIrrep {
    PermutationGroup group;
    Vector characters;

    XDIAG_API bool operator==(PermutationIrrep const &rhs) const;
    XDIAG_API bool operator!=(PermutationIrrep const &rhs) const;
  };

  XDIAG_API Representation() = default;
  XDIAG_API explicit Representation(PermutationGroup const &group);
  XDIAG_API Representation(PermutationGroup const &group,
                           Vector const &characters);
  XDIAG_API Representation(std::string type, int64_t charge);

  XDIAG_API bool operator==(Representation const &rhs) const;
  XDIAG_API bool operator!=(Representation const &rhs) const;

  // Discriminate the two variants without inspecting the type string.
  XDIAG_API bool is_permutation() const;
  XDIAG_API bool is_charge() const;

  XDIAG_API bool isreal() const;

  std::string type() const;
  PermutationGroup group() const;
  int64_t charge() const;
  Vector characters() const;

private:
  std::string type_; // always present; labels the variant
  std::variant<int64_t, PermutationIrrep> irrep_;
};

XDIAG_API bool isreal(Representation const &irrep);
XDIAG_API bool isapprox(Representation const &r1, Representation const &r2,
                        double tol = 1e-12);
XDIAG_API Representation multiply(Representation const &r1,
                                  Representation const &r2);
XDIAG_API Representation operator*(Representation const &r1,
                                   Representation const &r2);

XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   Representation const &irrep);
XDIAG_API std::string to_string(Representation const &irrep);
} // namespace xdiag
