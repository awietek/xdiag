#include "hashes.h"
#include <complex>
#include <hydra/random/hash_functions.h>

namespace hydra::random {

uint32_t hash(Permutation const &perm) {
  uint32_t h = 0;
  for (int i = 0; i < perm.size(); ++i) {
    h = hash_combine(h, hash_fnv1((uint32_t)perm[i]));
  }
  return h;
}

uint32_t hash(PermutationGroup const &group) {
  uint32_t h = 0;
  for (auto perm : group) {
    h = hash_combine(h, hash(perm));
  }
  return h;
}

uint32_t hash(Representation const &irrep) {
  uint32_t h = 0;
  for (auto sym : irrep.allowed_symmetries()) {
    h = hash_combine(h, hash_fnv1((uint32_t)sym));
  }
  for (auto ch : irrep.characters()) {
    double phase = std::arg(ch) + M_PI;
    double phase_scattered = std::log(
        std::abs(std::sin(1.7365529164217 * phase + 2.56381234623457)) +
        0.1234567);
    uint32_t angle = uint32_t(std::abs(phase_scattered * 1234567.0));
    h = hash_combine(h, hash_fnv1(angle));
  }
  return h;
}

template <typename bit_t> uint32_t hash(Spinhalf<bit_t> const &spinhalf) {
  uint32_t h =
      spinhalf.n_sites() == 0 ? 0 : hash_fnv1((uint32_t)spinhalf.n_sites());
  if (spinhalf.sz_conserved()) {
    h = hash_combine(h, hash_fnv1((uint32_t)spinhalf.n_up()));
  }
  if (spinhalf.symmetric()) {
    h = hash_combine(h, hash(spinhalf.permutation_group()));
    h = hash_combine(h, hash(spinhalf.irrep()));
  }
  return h;
}

template uint32_t hash(Spinhalf<uint16_t> const &);
template uint32_t hash(Spinhalf<uint32_t> const &);
template uint32_t hash(Spinhalf<uint64_t> const &);

template <typename bit_t> uint32_t hash(tJ<bit_t> const &tj) {
  uint32_t h = tj.n_sites() == 0 ? 0 : hash_fnv1((uint32_t)tj.n_sites());
  if (tj.charge_conserved() && tj.sz_conserved()) {
    h = hash_combine(h, hash_fnv1((uint32_t)tj.n_up()));
    h = hash_combine(h, hash_fnv1((uint32_t)tj.n_dn()));
  }
  if (tj.symmetric()) {
    h = hash_combine(h, hash(tj.permutation_group()));
    h = hash_combine(h, hash(tj.irrep()));
  }
  return h;
}
template uint32_t hash(tJ<uint16_t> const &);
template uint32_t hash(tJ<uint32_t> const &);
template uint32_t hash(tJ<uint64_t> const &);

template <typename bit_t> uint32_t hash(Electron<bit_t> const &electron) {
  uint32_t h =
      electron.n_sites() == 0 ? 0 : hash_fnv1((uint32_t)electron.n_sites());
  if (electron.charge_conserved() && electron.sz_conserved()) {
    h = hash_combine(h, hash_fnv1((uint32_t)electron.n_up()));
    h = hash_combine(h, hash_fnv1((uint32_t)electron.n_dn()));
  }
  if (electron.symmetric()) {
    h = hash_combine(h, hash(electron.permutation_group()));
    h = hash_combine(h, hash(electron.irrep()));
  }
  return h;
}

template uint32_t hash(Electron<uint16_t> const &);
template uint32_t hash(Electron<uint32_t> const &);
template uint32_t hash(Electron<uint64_t> const &);

} // namespace hydra::random
