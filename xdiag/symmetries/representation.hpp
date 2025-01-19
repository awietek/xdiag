#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <xdiag/common.hpp>

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/vector.hpp>

namespace xdiag {

class Representation {
public:
  XDIAG_API Representation() = default;
  XDIAG_API explicit Representation(PermutationGroup const &group);
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

  XDIAG_API int64_t size() const;
  XDIAG_API bool operator==(Representation const &rhs) const;
  XDIAG_API bool operator!=(Representation const &rhs) const;

  PermutationGroup const &group() const;
  Vector const &characters() const;
  bool isreal() const;

private:
  PermutationGroup group_;
  Vector characters_;
};

XDIAG_API bool isreal(Representation const &irrep);
XDIAG_API Representation multiply(Representation const &r1,
                                  Representation const &r2);
XDIAG_API Representation operator*(Representation const &r1,
                                   Representation const &r2);

XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   Representation const &irrep);
XDIAG_API std::string to_string(Representation const &irrep);
} // namespace xdiag
