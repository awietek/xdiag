#pragma once

#include <ostream>
#include <string>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

class Permutation {
public:
  XDIAG_API Permutation() = default;
  XDIAG_API Permutation(int64_t size); // creates identity permutation
  XDIAG_API Permutation(std::initializer_list<int64_t> list);
  XDIAG_API explicit Permutation(std::vector<int32_t> const &array);
  XDIAG_API explicit Permutation(std::vector<int64_t> const &array);
  XDIAG_API explicit Permutation(arma::Col<int64_t> const &array);
  XDIAG_API Permutation(int64_t *array, int64_t size);

  XDIAG_API int64_t size() const;
  XDIAG_API int64_t operator[](int64_t i) const;

  XDIAG_API Permutation &operator*=(Permutation const &rhs);
  XDIAG_API bool operator==(Permutation const &rhs) const;
  XDIAG_API bool operator!=(Permutation const &rhs) const;

  template <typename bit_t> bit_t apply(bit_t state) const;
  Permutation inverse() const;
  std::vector<int64_t> const &array() const;

private:
  std::vector<int64_t> array_;
};

XDIAG_API Permutation multiply(Permutation const &p1, Permutation const &p2);
XDIAG_API Permutation operator*(Permutation const &p1, Permutation const &p2);
XDIAG_API Permutation inverse(Permutation const &p);
XDIAG_API Permutation pow(Permutation const &p, int64_t power);

XDIAG_API std::ostream &operator<<(std::ostream &out, Permutation const &perm);
XDIAG_API std::string to_string(Permutation const &perm);

} // namespace xdiag
