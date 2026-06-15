// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>
#include <utility>

#include <xdiag/bits/nonzero.hpp>
#include <xdiag/bits/popcount.hpp>
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
// ("Cdagup"/"Cup" and "Cdagdn"/"Cdn") plus a cross sign. The decode is cold and
// lives in the .cpp (explicitly instantiated for the fermion/electron bit
// types); the per-state queries non_zero / action are hot and inline here.
template <typename bit_t> class CdagCString {
public:
  // Decode the `cdag_type` / `c_type` operators of `mono` into masks. Operators
  // of any other type are ignored (so a multi-sector monomial can be decoded
  // one sector at a time). Throws if the selected operators are not in
  // creation-left, strictly-ascending order.
  CdagCString(int64_t nsites, Monomial const &mono,
              std::string const &cdag_type, std::string const &c_type);

  // Sites whose occupation changes (number-operator sites cancel).
  bit_t flipmask() const;

  // Non-zero iff every C / number-operator site is occupied and every pure
  // Cdag site is empty.
  inline bool non_zero(bit_t spins) const {
    // pure Cdag sites = Cdag sites that are not also number-operator sites
    // (written with ^ and & only; bit_t need not provide operator~).
    bit_t mask_cdag_only = mask_cdag_ ^ (mask_cdag_ & mask_c_);
    return ((spins & mask_c_) == mask_c_) &&
           !bits::nonzero(spins & mask_cdag_only);
  }

  // {flipped configuration, negate}. Only meaningful when non_zero(spins). The
  // annihilation block acts first (each C seeing the initial occupancy below
  // it), then the creation block (each Cdag seeing spins ^ mask_c below it); the
  // Jordan-Wigner sign collapses to a single popcount. A single Cdag / C
  // reproduces the term_cdag / term_c kernels exactly.
  inline std::pair<bit_t, bool> action(bit_t spins) const {
    bit_t spins1 = spins ^ mask_c_;
    bit_t parity = (spins & signmask_c_) ^ (spins1 & signmask_cdag_);
    return {spins ^ (mask_c_ ^ mask_cdag_), (bool)(bits::popcount(parity) & 1)};
  }

private:
  bit_t mask_c_;        // sites carrying a C (must be occupied)
  bit_t mask_cdag_;     // sites carrying a Cdag (pure ones must be empty)
  bit_t signmask_c_;    // XOR of below(j_a) over the annihilation block
  bit_t signmask_cdag_; // XOR of below(i_b) over the creation block
};

} // namespace xdiag::matrices
