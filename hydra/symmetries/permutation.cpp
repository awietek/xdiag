#include "permutation.h"

#include <numeric>

#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

namespace hydra {

Permutation::Permutation(std::vector<int64_t> const &array)
    : size_(array.size()), array_(array) {

  // Check if permutation is valid
  for (int64_t i = 0; i < size_; ++i) {
    if (std::find(array_.begin(), array_.end(), i) == array_.end()) {
      Log.err("Error constructing Permutation: "
              "invalid permutation array");
    }
  }
}

Permutation::Permutation(std::initializer_list<int64_t> list)
    : Permutation(std::vector<int64_t>(list)) {}

Permutation::Permutation(io::FileTomlHandler && hdl)
    : Permutation(hdl.as<Permutation>()) {}

template <typename bit_t> bit_t Permutation::apply(bit_t state) const {
  bit_t tstate = 0;
  for (int64_t site = 0; site < size_; ++site) {
    tstate |= ((state >> site) & 1) << array_[site];
  }
  return tstate;
}

template uint16_t Permutation::apply<uint16_t>(uint16_t state) const;
template uint32_t Permutation::apply<uint32_t>(uint32_t state) const;
template uint64_t Permutation::apply<uint64_t>(uint64_t state) const;

Permutation Permutation::inverse() const {
  std::vector<int64_t> perm_inv(size_, 0);
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

int64_t Permutation::size() const { return size_; }
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

Permutation operator*(Permutation const &p1, Permutation const &p2) {
  if (p1.size() != p2.size()) {
    Log.err("Error multiplying Permutation: the two permutations do not have "
            "the same number of sites. p1.size()={}, p2.size()={}",
            p1.size(), p2.size());
  }
  int64_t size = p1.size();
  std::vector<int64_t> array(size, 0);
  for (int64_t i = 0; i < size; ++i) {
    array[i] = p1[p2[i]];
  }
  return Permutation(array);
}

Permutation inverse(Permutation const &p) { return p.inverse(); }
Permutation shuffle(Permutation const &p) { return p.shuffle(); }

} // namespace hydra
