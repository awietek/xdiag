#include "permutation.hpp"

#include <numeric>

#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

static void check_valid_permutation(std::vector<int64_t> const &array) try {
  for (int64_t i = 0; i < array.size(); ++i) {
    if (std::find(array.begin(), array.end(), i) == array.end()) {
      XDIAG_THROW("Invalid permutation array. A Permutation is valid if every "
                  "index from 0 to size-1 is present exactly once.");
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation::Permutation(int64_t size) : array_(size, 0) {
  std::iota(array_.begin(), array_.end(), 0);
}

Permutation::Permutation(std::initializer_list<int64_t> list) try
    : Permutation(std::vector<int64_t>(list)) {
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation::Permutation(std::vector<int32_t> const &array) try
    : array_(array.size()) {
  int64_t idx = 0;
  for (int32_t p : array) {
    array_[idx] = (int64_t)p;
    ++idx;
  }
  check_valid_permutation(array_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation::Permutation(std::vector<int64_t> const &array) try
    : array_(array) {
  check_valid_permutation(array_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation::Permutation(arma::Col<int64_t> const &array) try
    : array_(array.memptr(), array.memptr() + array.size()) {
  check_valid_permutation(array_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation::Permutation(int64_t *ptr, int64_t size) try
    : array_(ptr, ptr + size) {
  check_valid_permutation(array_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename bit_t> bit_t Permutation::apply(bit_t state) const {
  bit_t tstate = 0;
  for (int64_t site = 0; site < array_.size(); ++site) {
    tstate |= ((state >> site) & 1) << array_[site];
  }
  return tstate;
}

template uint16_t Permutation::apply<uint16_t>(uint16_t state) const;
template uint32_t Permutation::apply<uint32_t>(uint32_t state) const;
template uint64_t Permutation::apply<uint64_t>(uint64_t state) const;

Permutation Permutation::inverse() const try {
  std::vector<int64_t> perm_inv(array_.size(), 0);
  int64_t idx = 0;
  for (auto p : array_) {
    perm_inv[p] = idx;
    ++idx;
  }
  return Permutation(perm_inv);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation &Permutation::operator*=(Permutation const &rhs) try {
  if (array_.size() != rhs.size()) {
    XDIAG_THROW(
        fmt::format("The two permutations do not have "
                    "the same number of sites. p1.size()={}, p2.size()={}",
                    array_.size(), rhs.size()));
  }
  int64_t size = array_.size();
  std::vector<int64_t> array(size, 0);
  for (int64_t i = 0; i < size; ++i) {
    array[i] = array_[rhs[i]];
  }
  std::swap(array_, array);
  return *this;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t Permutation::size() const { return array_.size(); }
int64_t Permutation::operator[](int64_t i) const { return array_[i]; }
bool Permutation::operator==(Permutation const &rhs) const {
  return rhs.array_ == array_;
}
bool Permutation::operator!=(Permutation const &rhs) const {
  return !operator==(rhs);
}

std::vector<int64_t> const &Permutation::array() const { return array_; }

int64_t size(Permutation const &p) { return p.size(); }
Permutation multiply(Permutation const &p1, Permutation const &p2) try {
  Permutation p = p1;
  p *= p2;
  return p;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation operator*(Permutation const &p1, Permutation const &p2) try {
  return multiply(p1, p2);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Permutation inverse(Permutation const &p) { return p.inverse(); }
XDIAG_API Permutation pow(Permutation const &p, int64_t power) try {
  Permutation pp(p.size());
  if (power >= 0) {
    for (int i = 0; i < power; ++i) {
      pp *= p;
    }
  } else {
    Permutation pi = inverse(p);
    for (int i = 0; i < -power; ++i) {
      pp *= pi;
    }
  }
  return pp;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::ostream &operator<<(std::ostream &out, Permutation const &p) {
  for (int64_t i = 0; i < p.size(); ++i) {
    out << p[i] << " ";
  }
  return out;
}
std::string to_string(Permutation const &perm) {
  return to_string_generic(perm);
}

} // namespace xdiag
