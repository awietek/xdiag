#pragma once

#include <vector>

namespace hydra {

class Permutation {
public:
  Permutation() = default;
  explicit Permutation(std::vector<int> const &perm);

  template <typename bit_t> bit_t apply(bit_t state) const;
  Permutation inverse() const;
  Permutation shuffle() const;

  inline int n_sites() const { return n_sites_; }
  inline int operator[](int i) const { return permutation_[i]; }
  inline bool operator==(Permutation const &rhs) const {
    return rhs.permutation_ == permutation_;
  }
  inline bool operator!=(Permutation const &rhs) const {
    return !operator==(rhs);
  }

private:
  int n_sites_;
  std::vector<int> permutation_;
};

Permutation IdentityPermutation(int n_sites);
Permutation operator*(Permutation const &p1, Permutation const &p2);
Permutation Inverse(Permutation const &p);
Permutation Shuffle(Permutation const &p);

} // namespace hydra
