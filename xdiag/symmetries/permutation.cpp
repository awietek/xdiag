#include "permutation.hpp"

#include <numeric>

#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

void check_valid_permutation(std::vector<int64_t> const &array) {
  for (int64_t i = 0; i < array.size(); ++i) {
    if (std::find(array.begin(), array.end(), i) == array.end()) {
      XDIAG_THROW("Error constructing Permutation: "
                  "invalid permutation array");
    }
  }
}

Permutation::Permutation(std::vector<int32_t> const &array)
    : array_(array.size()) {
  int64_t idx = 0;
  for (int32_t p : array) {
    array_[idx] = p;
    ++idx;
  }
  check_valid_permutation(array_);
}

Permutation::Permutation(std::vector<int64_t> const &array) : array_(array) {
  check_valid_permutation(array_);
}

Permutation::Permutation(io::FileTomlHandler &&hdl)
    : Permutation(hdl.as<Permutation>()) {}

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

Permutation Permutation::inverse() const {
  std::vector<int64_t> perm_inv(array_.size(), 0);
  int64_t idx = 0;
  for (auto p : array_) {
    perm_inv[p] = idx;
    ++idx;
  }
  return Permutation(perm_inv);
}

Permutation Permutation::shuffle() const {
  std::vector<int64_t> ps = array_;
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(ps.begin(), ps.end(), g);
  return Permutation(ps);
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

Permutation identity_permutation(int64_t size) {
  std::vector<int64_t> array(size, 0);
  std::iota(array.begin(), array.end(), 0);
  return Permutation(array);
}

Permutation operator*(Permutation const &p1, Permutation const &p2) try {
  if (p1.size() != p2.size()) {
    XDIAG_THROW(fmt::format(
        "Error multiplying Permutation: the two permutations do not have "
        "the same number of sites. p1.size()={}, p2.size()={}",
        p1.size(), p2.size()));
  }
  int64_t size = p1.size();
  std::vector<int64_t> array(size, 0);
  for (int64_t i = 0; i < size; ++i) {
    array[i] = p1[p2[i]];
  }
  return Permutation(array);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return Permutation();
}

Permutation inverse(Permutation const &p) { return p.inverse(); }
Permutation shuffle(Permutation const &p) { return p.shuffle(); }

} // namespace xdiag
