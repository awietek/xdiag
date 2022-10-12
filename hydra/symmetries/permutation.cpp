#include "permutation.h"

#include <numeric>

#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

namespace hydra {

Permutation::Permutation(std::vector<int> const &perm)
    : n_sites_(perm.size()), permutation_(perm) {

  // Check if permutation is valid
  for (int i = 0; i < n_sites_; ++i) {
    if (std::find(permutation_.begin(), permutation_.end(), i) ==
        permutation_.end()) {
      Log.err("Error constructing Permutation: "
              "invalid permutation array");
    }
  }
}

Permutation::Permutation(std::initializer_list<int> perm)
    : Permutation(std::vector<int>(perm)) {}

template <typename bit_t> bit_t Permutation::apply(bit_t state) const {
  bit_t tstate = 0;
  for (int site = 0; site < n_sites_; ++site) {
    tstate |= ((state >> site) & 1) << permutation_[site];
  }
  return tstate;
}

template uint16_t Permutation::apply<uint16_t>(uint16_t state) const;
template uint32_t Permutation::apply<uint32_t>(uint32_t state) const;
template uint64_t Permutation::apply<uint64_t>(uint64_t state) const;

Permutation Permutation::inverse() const {
  std::vector<int> perm_inv(n_sites_, 0);
  int idx = 0;
  for (auto p : permutation_) {
    perm_inv[p] = idx;
    ++idx;
  }
  return Permutation(perm_inv);
}

Permutation Permutation::shuffle() const {
  std::vector<int> ps = permutation_;
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(ps.begin(), ps.end(), g);
  return Permutation(ps);
}

Permutation identity_permutation(int n_sites) {
  std::vector<int> perm(n_sites, 0);
  std::iota(perm.begin(), perm.end(), 0);
  return Permutation(perm);
}

Permutation operator*(Permutation const &p1, Permutation const &p2) {
  if (p1.n_sites() != p2.n_sites()) {
    Log.err("Error multiplying Permutation: the two permutations do not have "
            "the same number of sites. p1.n_sites()={}, p2.n_sites()={}",
            p1.n_sites(), p2.n_sites());
  }
  int n_sites = p1.n_sites();
  std::vector<int> perm(n_sites, 0);
  for (int i = 0; i < n_sites; ++i) {
    perm[i] = p1[p2[i]];
  }
  return Permutation(perm);
}

Permutation inverse(Permutation const &p) { return p.inverse(); }
Permutation shuffle(Permutation const &p) { return p.shuffle(); }

} // namespace hydra
