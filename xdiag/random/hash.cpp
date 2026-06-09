// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hash.hpp"

#include <cmath>
#include <complex>
#include <variant>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/numbers.hpp>
#include <xdiag/random/hash_functions.hpp>

namespace xdiag::random {

uint64_t hash(Permutation const &perm) {
  uint64_t h = 0;
  for (int64_t i = 0; i < perm.size(); ++i) {
    h = hash_combine(h, hash_fnv1((uint64_t)perm[i]));
  }
  return h;
}

uint64_t hash(PermutationGroup const &group) {
  uint64_t h = 0;
  for (int64_t i = 0; i < group.size(); ++i) {
    h = hash_combine(h, hash(group[i]));
  }
  return h;
}

// FNV-1 hash of a string, folding in its bytes. Used so Representations of
// different type (e.g. "number" vs "nup") hash differently even at equal charge.
static uint64_t hash_string(std::string const &s) {
  uint64_t h = 0;
  for (char c : s) {
    h = hash_combine(h, hash_fnv1((uint64_t)(unsigned char)c));
  }
  return h;
}

uint64_t hash(Representation const &irrep) {
  uint64_t h = hash_string(irrep.type());
  if (irrep.is_charge()) { // U(1) charge irrep (e.g. "number", "nup")
    h = hash_combine(h, hash_fnv1((uint64_t)irrep.charge()));
  } else { // permutation irrep: fold in characters and the group
    for (std::complex<double> ch : irrep.characters().as<arma::cx_vec>()) {
      double phase = std::arg(ch) + math::pi;
      double phase_scattered = std::log(
          std::abs(std::sin(1.7365529164217 * phase + 2.56381234623457)) +
          0.1234567);
      uint64_t angle = uint64_t(std::abs(phase_scattered * 1234567.0));
      h = hash_combine(h, hash_fnv1(angle));
    }
    h = hash_combine(h, hash(irrep.group()));
  }
  return h;
}

uint64_t hash(Block const &block) {
  return std::visit([](auto &&block) { return hash(block); }, block);
}

uint64_t hash(Spinhalf const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  // Fold in every irrep the block carries (U(1) charges such as "nup", and a
  // "SitePermutation" irrep). XOR keeps this order-independent, matching the
  // order-independent equality of RepresentationSet.
  uint64_t irreps_hash = 0;
  for (Representation const &irrep : block.irreps()) {
    irreps_hash ^= hash(irrep);
  }
  return hash_combine(h, irreps_hash);
}

uint64_t hash(Boson const &block) {
  uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
  h = hash_combine(h, hash_fnv1((uint64_t)block.d()));
  // Fold in every irrep the block carries (U(1) charges such as "number", and a
  // "SitePermutation" irrep). XOR keeps this order-independent, matching the
  // order-independent equality of RepresentationSet.
  uint64_t irreps_hash = 0;
  for (Representation const &irrep : block.irreps()) {
    irreps_hash ^= hash(irrep);
  }
  return hash_combine(h, irreps_hash);
}

// uint64_t hash(tJ const &block) {
//   uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
//   if (block.nup()) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
//   }
//   if (block.ndn()) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
//   }
//   if (block.irrep()) {
//     h = hash_combine(h, hash(*block.irrep()));
//   }
//   return h;
// }

// uint64_t hash(Electron const &block) {
//   uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
//   if (block.nup()) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
//   }
//   if (block.ndn()) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
//   }
//   if (block.irrep()) {
//     h = hash_combine(h, hash(*block.irrep()));
//   }
//   return h;
// }

// #ifdef XDIAG_USE_MPI

// uint64_t hash(SpinhalfDistributed const &block) {
//   uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
//   if (block.nup() != undefined) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
//   }
//   int mpi_rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//   h = hash_combine(h, hash_fnv1((uint64_t)mpi_rank));
//   return h;
// }

// uint64_t hash(tJDistributed const &block) {
//   uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
//   if (block.nup() != undefined) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
//   }
//   if (block.ndn() != undefined) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
//   }

//   int mpi_rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//   h = hash_combine(h, hash_fnv1((uint64_t)mpi_rank));
//   return h;
// }

// uint64_t hash(ElectronDistributed const &block) {
//   uint64_t h = block.nsites() == 0 ? 0 : hash_fnv1((uint64_t)block.nsites());
//   if (block.nup() != undefined) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.nup()));
//   }
//   if (block.ndn() != undefined) {
//     h = hash_combine(h, hash_fnv1((uint64_t)*block.ndn()));
//   }

//   int mpi_rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//   h = hash_combine(h, hash_fnv1((uint64_t)mpi_rank));
//   return hash_combine(h, 1234);
// }
// #endif

} // namespace xdiag::random
