// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "permutation.hpp"

#include <numeric>

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>
#include <xdiag/utils/xdiag_offset.hpp>

namespace xdiag {

static void check_valid_permutation(std::vector<int64_t> const &array) try {
  for (int64_t i = 0; i < array.size(); ++i) {
    if (std::find(array.begin(), array.end(), i) == array.end()) {
      XDIAG_THROW("Invalid permutation array. A Permutation is valid if every "
                  "index from 0 to size-1 is present exactly once.");
    }
  }
}
XDIAG_CATCH

Permutation::Permutation(int64_t size) : array_(size, 0) {
  std::iota(array_.begin(), array_.end(), 0);
}

Permutation::Permutation(std::initializer_list<int64_t> list) try
    : Permutation(std::vector<int64_t>(list)) {}
XDIAG_CATCH

Permutation::Permutation(std::vector<int32_t> const &array) try
    : array_(array.size()) {
  int64_t idx = 0;
  for (int32_t p : array) {
    array_[idx] = (int64_t)p;
    ++idx;
  }
  check_valid_permutation(array_);
}
XDIAG_CATCH

Permutation::Permutation(std::vector<int64_t> const &array) try
    : array_(array) {
  check_valid_permutation(array_);
}
XDIAG_CATCH

Permutation::Permutation(arma::Col<int64_t> const &array) try
    : array_(array.memptr(), array.memptr() + array.size()) {
  check_valid_permutation(array_);
}
XDIAG_CATCH

Permutation::Permutation(int64_t const *ptr, int64_t size) try
    : array_(ptr, ptr + size) {
  check_valid_permutation(array_);
}
XDIAG_CATCH

Permutation Permutation::inv() const try {
  std::vector<int64_t> perm_inv(array_.size(), 0);
  int64_t idx = 0;
  for (auto p : array_) {
    perm_inv[p] = idx;
    ++idx;
  }
  return Permutation(perm_inv);
}
XDIAG_CATCH

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
}
XDIAG_CATCH

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
}
XDIAG_CATCH

Permutation operator*(Permutation const &p1, Permutation const &p2) try {
  return multiply(p1, p2);
}
XDIAG_CATCH

Permutation inv(Permutation const &p) { return p.inv(); }
Permutation pow(Permutation const &p, int64_t power) try {
  Permutation pp(p.size());
  if (power >= 0) {
    for (int i = 0; i < power; ++i) {
      pp *= p;
    }
  } else {
    Permutation pi = inv(p);
    for (int i = 0; i < -power; ++i) {
      pp *= pi;
    }
  }
  return pp;
}
XDIAG_CATCH

std::ostream &operator<<(std::ostream &out, Permutation const &p) {
  for (int64_t i = 0; i < p.size(); ++i) {
    out << p[i] + XDIAG_OFFSET << " ";
  }
  return out;
}
std::string to_string(Permutation const &perm) {
  return utils::to_string_generic(perm);
}

} // namespace xdiag
