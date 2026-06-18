// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>
#include <utility>

#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/operators/monomial.hpp>

namespace xdiag::matrices {

// State-independent encoding of a normal-ordered single-sector fermion operator
// string acting on one bit configuration. Every operator of `mono` must be of
// type `cdag_type` or `c_type` (anything else throws), and they are assumed to
// be in the textbook creation-left normal order produced by the normal-order
// rules:
//
//   Cdag{i_1} ... Cdag{i_p}  C{j_1} ... C{j_q} ,   i_1<...<i_p ,  j_1<...<j_q
//
// Same-site nilpotency / CAR is already resolved, so a site carries at most one
// C, one Cdag, or the pair Cdag*C (a number operator).
//
// This is the shared core of the fermion and electron Cdag/C string kernels:
// the fermion block builds one string ("Cdag", "C"); the electron block splits
// its monomial into an up and a dn sub-monomial and builds one string each
// ("Cdagup"/"Cup" and "Cdagdn"/"Cdn") plus a cross sign. 
template <typename bit_t> class CdagCString {
public:
  // Decode the `cdag_type` / `c_type` operators of `mono` into masks. Every
  // operator of `mono` must be of one of those two types (a multi-sector
  // monomial is split by the caller first). Throws on any other type, or if the
  // operators are not in creation-left, strictly-ascending order.
  CdagCString(int64_t nsites, Monomial const &mono,
              std::string const &cdag_type, std::string const &c_type);

  // Non-zero iff every C / number-operator site is occupied and every pure
  // Cdag site is empty.
  inline bool non_zero(bit_t spins) const {
    return ((spins & mask_c_) == mask_c_) &&
           bits::iszero(spins & mask_cdag_only_);
  }

  // {flipped configuration, negate}. Only meaningful when non_zero(spins). The
  // annihilation block acts first (each C seeing the initial occupancy below
  // it), then the creation block (each Cdag seeing spins ^ mask_c below it); the
  // Jordan-Wigner sign collapses to a single popcount. A single Cdag / C
  // reproduces the term_cdag / term_c kernels exactly.
  inline std::pair<bit_t, bool> action(bit_t spins) const {
    bit_t spins1 = spins ^ mask_c_;
    bit_t parity = (spins & signmask_c_) ^ (spins1 & signmask_cdag_);
    return {spins ^ flipmask_, (bool)(bits::popcount(parity) & 1)};
  }

private:
  bit_t mask_c_;         // sites carrying a C (must be occupied)
  bit_t mask_cdag_only_; // pure Cdag sites, i.e. not number-operator sites
                         // (must be empty)
  bit_t flipmask_;       // sites whose occupation changes (mask_c ^ mask_cdag)
  bit_t signmask_c_;     // XOR of below(j_a) over the annihilation block
  bit_t signmask_cdag_;  // XOR of below(i_b) over the creation block
};

// Electron / tJ string kernels split a two-sector monomial into an up and a dn
// sub-string and evaluate them independently, which silently reorders the
// operators into the all-ups-then-all-dns form. That reordering carries a
// Jordan-Wigner sign whenever a dn operator originally sat to the LEFT of an up
// operator (e.g. S- = Cdagdn_i Cup_i, or the interleaved strings the tJ
// normal-order rule emits because it cannot separate the projected S+/S- pairs).
// Returns true iff that partition sign is negative: one factor of -1 per
// (dn-operator before up-operator) pair. It is +1 for an already separated
// monomial (the electron normal order), so folding it into the coefficient
// leaves the electron block unchanged while making the kernels correct for any
// operator order. `cdag_dn` / `c_dn` name the dn-sector operator types.
inline bool cdagc_sector_partition_neg(Monomial const &mono,
                                       std::string const &cdag_dn,
                                       std::string const &c_dn) {
  int64_t dn_seen = 0;
  int64_t inversions = 0;
  for (int64_t k = 0; k < mono.size(); ++k) {
    std::string const &t = mono[k].type();
    if ((t == cdag_dn) || (t == c_dn)) {
      ++dn_seen; // a dn operator: later up operators must hop past it
    } else {
      inversions += dn_seen;
    }
  }
  return (inversions & 1) != 0;
}

} // namespace xdiag::matrices
