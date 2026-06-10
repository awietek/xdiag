// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cassert>
#include <string>
#include <utility>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/nonzero.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::matrices::fermion {

// Applies a normal-ordered string of elementary fermion operators to a fermion
// basis. The string is assumed to come out of the fermion algebra in textbook
// normal order:
//
//   Cdag{i_1} ... Cdag{i_p}  C{j_1} ... C{j_q} ,   i_1<...<i_p ,  j_1<...<j_q
//
// i.e. all creation operators (ascending site) to the left of all annihilation
// operators (ascending site). Same-site nilpotency and the CAR have already
// been resolved, so a site carries at most one C, one Cdag, or the pair
// Cdag*C (the number operator, which is normal ordered).
//
// Decoding into state-independent masks:
//
//   mask_c        sites carrying a C    (must be occupied)
//   mask_cdag     sites carrying a Cdag (pure ones must be empty)
//   flipmask    = mask_c ^ mask_cdag      sites whose occupation changes
//                                         (number-operator sites cancel)
//   signmask_c    = XOR_a below(j_a)      JW mask of the annihilation block
//   signmask_cdag = XOR_b below(i_b)      JW mask of the creation block
//   below(s)      = bits with index < s
//
// Sign. The annihilation block acts first (rightmost), each C seeing the
// *initial* occupancy below it (ascending block, no lower operator to its
// right). It then clears every mask_c site, leaving spins1 = spins ^ mask_c.
// The creation block acts next, each Cdag seeing the occupancy of spins1 below
// it. Hence
//
//   sign = (-1)^[ popcount(spins  & signmask_c)
//               + popcount(spins1 & signmask_cdag) ] ,   spins1 = spins ^ mask_c
//
// collapsed to a single popcount via parity linearity. This uses the same
// "occupied sites with index < s" Jordan-Wigner convention as term_c /
// term_cdag, so a single Cdag / C reproduces those kernels exactly.
template <typename coeff_t, class basis_t, class fill_f>
void term_cdagc_string(Coeff const &c, Monomial const &mono,
                       basis_t const &basis_in, basis_t const &basis_out,
                       fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();

  // Decode the string once into state-independent masks.
  bit_t mask_c = bits::zero<bit_t>(nsites);
  bit_t mask_cdag = bits::zero<bit_t>(nsites);
  bit_t signmask_c = bits::zero<bit_t>(nsites);
  bit_t signmask_cdag = bits::zero<bit_t>(nsites);
  bool in_c_block = false;       // flips true once the C block starts
  int64_t prev_cdag_site = -1;   // for the ascending-order assert
  int64_t prev_c_site = -1;
  for (int64_t k = 0; k < mono.size(); ++k) {
    Op const &op = mono[k];
    int64_t s = op[0];
    std::string type = op.type();
    bit_t below = bits::bitmask<bit_t>(nsites, s); // bits with index < s

    if (type == "Cdag") {
      // Creation-left form: no Cdag may follow a C, and Cdag sites within the
      // block must be strictly ascending. Required for the sign to be correct.
      assert(!in_c_block &&
             "term_cdagc_string expects creation-left normal order: a Cdag "
             "follows a C. Normal-order through the fermion algebra first.");
      assert((s > prev_cdag_site) &&
             "term_cdagc_string expects the Cdag block strictly ascending.");
      prev_cdag_site = s;
      bits::set(mask_cdag, s);
      signmask_cdag ^= below;
    } else if (type == "C") {
      in_c_block = true;
      assert((s > prev_c_site) &&
             "term_cdagc_string expects the C block strictly ascending.");
      prev_c_site = s;
      bits::set(mask_c, s);
      signmask_c ^= below;
    } else {
      XDIAG_THROW(fmt::format(
          "Invalid Op type \"{}\" encountered in fermion Cdag/C string; only "
          "\"Cdag\" and \"C\" are allowed.",
          type));
    }
  }

  bit_t flipmask = mask_c ^ mask_cdag;
  // pure creation sites = Cdag sites that are not also number-operator sites.
  // Written with ^ and & only (avoids operator~, which not every bit_t has).
  bit_t mask_cdag_only = mask_cdag ^ (mask_cdag & mask_c);

  // Non-zero iff every C/number site is occupied and every pure Cdag site is
  // empty.
  auto non_zero_term = [=](bit_t spins) -> bool {
    return ((spins & mask_c) == mask_c) &&
           !bits::nonzero(spins & mask_cdag_only);
  };
  auto term_action = [=](bit_t spins) -> std::pair<bit_t, coeff_t> {
    // spins1: state after the annihilation block has cleared all mask_c sites
    // (valid here, since non_zero_term guarantees they are all occupied).
    bit_t spins1 = spins ^ mask_c;
    bit_t parity = (spins & signmask_c) ^ (spins1 & signmask_cdag);
    coeff_t sign = (bits::popcount(parity) & 1) ? -cf : cf;
    return {spins ^ flipmask, sign};
  };
  term_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::fermion
