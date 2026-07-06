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

class XDIAG_API Representation {
public:
  // A 1-D irrep of a finite group of site permutations: one character per
  // group element.
  struct PermutationIrrep {
    PermutationGroup group;
    Vector characters;

    bool operator==(PermutationIrrep const &rhs) const;
    bool operator!=(PermutationIrrep const &rhs) const;
  };

  Representation() = default;
  explicit Representation(PermutationGroup const &group);
  Representation(PermutationGroup const &group, Vector const &characters);
  Representation(std::string type, int64_t charge);

  bool operator==(Representation const &rhs) const;
  bool operator!=(Representation const &rhs) const;

  // Discriminate the two variants without inspecting the type string.
  bool is_permutation() const;
  bool is_charge() const;

  bool isreal() const;

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
