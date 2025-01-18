#pragma once

#include <vector>

#include <xdiag/common.hpp>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/vector.hpp>

namespace xdiag {

class Representation {
public:
  XDIAG_API Representation() = default;

  template <typename T>
  XDIAG_API Representation(PermutationGroup const &group,
                           std::vector<T> const &characters);
  template <typename T>
  XDIAG_API Representation(PermutationGroup const &group,
                           arma::Col<T> const &characters);
  template <typename T>
  XDIAG_API Representation(PermutationGroup const &group, T *characters,
                           int64_t n_characters);

  Representation(PermutationGroup const &group, Vector const &characters);

  PermutationGroup const &group() const;
  Vector const &characters() const;

  int64_t size() const;
  bool isreal() const;
  bool operator==(Representation const &rhs) const;
  bool operator!=(Representation const &rhs) const;

private:
  PermutationGroup group_;
  Vector characters_;
};

bool isreal(Representation const &irrep);

Representation trivial_representation(int64_t size);
Representation trivial_representation(PermutationGroup const &group);

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep);

Representation multiply(Representation const &r1, Representation const &r2);
Representation operator*(Representation const &r1, Representation const &r2);

std::ostream &operator<<(std::ostream &out, Representation const &irrep);
std::string to_string(Representation const &irrep);
} // namespace xdiag
