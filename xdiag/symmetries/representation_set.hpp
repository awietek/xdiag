// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <initializer_list>
#include <optional>
#include <string>
#include <vector>

#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

// A collection of Representations, each identified by its type string. Although
// backed by a std::vector, the type acts as a key: a RepresentationSet holds at
// most one Representation per type and compares equal independent of the order
// in which the Representations were given.
class RepresentationSet {
public:
  RepresentationSet() = default;
  RepresentationSet(std::vector<Representation> const &representations);
  RepresentationSet(std::initializer_list<Representation> representations);

  bool has_type(std::string type) const;
  // Return the requested quantity for the Representation of the given type, or
  // std::nullopt if no such Representation is present (or it does not carry the
  // requested quantity).
  std::optional<PermutationGroup> group(std::string type) const;
  std::optional<int64_t> charge(std::string type) const;
  std::optional<Vector> characters(std::string type) const;

  bool operator==(RepresentationSet const &rhs) const;
  bool operator!=(RepresentationSet const &rhs) const;

  bool isreal() const;
  int64_t size() const;

  bool isapprox(RepresentationSet const &rhs, double tol = 1e-12) const;
  RepresentationSet multiply(RepresentationSet const &rhs) const;

  // Iteration over the contained Representations (in insertion order).
  using const_iterator = std::vector<Representation>::const_iterator;
  const_iterator begin() const;
  const_iterator end() const;

private:
  std::vector<Representation> representations_;
};

bool isreal(RepresentationSet const &set);
bool isapprox(RepresentationSet const &s1, RepresentationSet const &s2,
              double tol = 1e-12);

RepresentationSet multiply(RepresentationSet const &s1,
                           RepresentationSet const &s2);
RepresentationSet operator*(RepresentationSet const &s1,
                            RepresentationSet const &s2);

} // namespace xdiag
