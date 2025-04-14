// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hash.hpp"

#include <complex>
#include <variant>

#include <xdiag/random/hash_functions.hpp>

namespace xdiag::random {

uint64_t hash(Permutation const &perm) {
  uint64_t h = 0;
  for (int i = 0; i < perm.size(); ++i) {
    h = hash_combine(h, hash_fnv1((uint64_t)perm[i]));
  }
  return h;
}

uint64_t hash(PermutationGroup const &group) {
  uint64_t h = 0;
  for (auto perm : group) {
    h = hash_combine(h, hash(perm));
  }
  return h;
}

uint64_t hash(Representation const &irrep) {
  uint64_t h = 0;
  for (auto ch : irrep.characters().as<arma::cx_vec>()) {
    double phase = std::arg(ch) + XDIAG_PI;
    double phase_scattered = std::log(
        std::abs(std::sin(1.7365529164217 * phase + 2.56381234623457)) +
        0.1234567);
    uint64_t angle = uint64_t(std::abs(phase_scattered * 1234567.0));
    h = hash_combine(h, hash_fnv1(angle));
  }
  return hash_combine(h, hash(irrep.group()));
}

uint64_t hash(Block const &block) {
  return std::visit([](auto &&block) { return hash(block); }, block);
}

uint64_t hash(Spinhalf const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  if (block.nup()) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
  }
  if (block.irrep()) {
    h = hash_combine(h, hash(*block.irrep()));
  }
  return h;
}

uint64_t hash(tJ const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  if (block.nup()) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
  }
  if (block.ndn()) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
  }
  if (block.irrep()) {
    h = hash_combine(h, hash(*block.irrep()));
  }
  return h;
}

uint64_t hash(Electron const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  if (block.nup()) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
  }
  if (block.ndn()) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
  }
  if (block.irrep()) {
    h = hash_combine(h, hash(*block.irrep()));
  }
  return h;
}

#ifdef XDIAG_USE_MPI

uint64_t hash(SpinhalfDistributed const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  if (block.nup() != undefined) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
  }
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  h = hash_combine(h, hash_fnv1((uint64_t)mpi_rank));
  return h;
}

uint64_t hash(tJDistributed const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  if (block.nup() != undefined) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
  }
  if (block.ndn() != undefined) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
  }

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  h = hash_combine(h, hash_fnv1((uint64_t)mpi_rank));
  return h;
}

uint64_t hash(ElectronDistributed const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  if (block.nup() != undefined) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
  }
  if (block.ndn() != undefined) {
    h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
  }

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  h = hash_combine(h, hash_fnv1((uint64_t)mpi_rank));
  return hash_combine(h, 1234);
}
#endif

} // namespace xdiag::random
