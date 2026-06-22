// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/kernels/fill_functions.hpp>
#include <xdiag/kernels/terms/cdagc_string.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::kernels::tj {

// General product of elementary operators (Cdagup/Cup/Cdagdn/Cdn) on the tJ
// basis -- the cold fallback / correctness oracle for the fast named kernels.
//
// The tJ space is the electron space projected to no double occupancy, so the
// matrix elements (and their Jordan-Wigner signs) are exactly the electron ones
// restricted to that subspace. We therefore reuse the electron machinery
// verbatim: split the monomial into an up sub-string and a dn sub-string, apply
// CdagCString to each (full configurations), and add the cross sign (every dn
// operator sits behind the whole up string -> (-1)^(q_dn * Nup_in)). The only
// tJ-specific additions are: decompress the dn before applying dn_str, drop any
// output that is doubly occupied (ups_mid & dns_mid != 0 -- this also correctly
// handles strings that move an up and a dn through the same site, e.g.
// Exchange), and recompress for the index.
//
// Unlike the named kernels this materialises the dn configuration (an O(nsites)
// deposit/extract per non-zero), so it is deliberately confined to the cold path
// once Hop*/Exchange and the raise/lower ops are routed to their O(1) kernels.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagc_string(Coeff const &c, Monomial const &mono,
                       basis::BasistJ<enumeration_t> const &basis_in,
                       basis::BasistJ<enumeration_t> const &basis_out,
                       fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  // The split below evaluates the up and dn sub-strings as if the monomial were
  // in all-ups-then-all-dns order; fold in the Jordan-Wigner sign of partitioning
  // the (possibly interleaved) input into that order.
  if (cdagc_sector_partition_neg(mono, "Cdagdn", "Cdn")) {
    cf = -cf;
  }

  // Split the monomial into its up and dn sub-strings (validating the types).
  std::vector<Op> up_ops, dn_ops;
  for (int64_t k = 0; k < mono.size(); ++k) {
    Op const &op = mono[k];
    std::string type = op.type();
    if ((type == "Cdagup") || (type == "Cup")) {
      up_ops.push_back(op);
    } else if ((type == "Cdagdn") || (type == "Cdn")) {
      dn_ops.push_back(op);
    } else {
      XDIAG_THROW("tJ Cdag/C string: unexpected operator type \"" + type +
                  "\"; only Cdagup, Cup, Cdagdn, Cdn are allowed");
    }
  }
  CdagCString<bit_t> up_str(nsites, Monomial(up_ops), "Cdagup", "Cup");
  CdagCString<bit_t> dn_str(nsites, Monomial(dn_ops), "Cdagdn", "Cdn");
  bool cross_when_odd_nup = (dn_ops.size() & 1);

  auto const &basis_up = basis_in.basis_up();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up, num_thread, omp_get_num_threads());
#else
    auto [begin_up, end_up, idx_up] = utils::thread_range(basis_up, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      if (!up_str.non_zero(ups)) {
        continue;
      }
      // up sub-string action + cross sign, hoisted out of the dn loop
      auto [ups_mid, up_total] = up_str.action(ups);
      if (cross_when_odd_nup && (bits::popcount(ups) & 1)) {
        up_total = !up_total;
      }
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(basis_out.index_up(ups_mid));

      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        bit_t dns = basis::tj_decompress_dns(ups, dnc, nsites);
        if (dn_str.non_zero(dns)) {
          auto [dns_mid, dn_fermi] = dn_str.action(dns);
          if (bits::iszero(ups_mid & dns_mid)) { // no double occupancy
            bit_t dnc_out = basis::tj_compress_dns(ups_mid, dns_mid, nsites);
            int64_t idx_out = base_out + basis_out.index_dncs(ups_mid, dnc_out);
            bool neg = up_total ^ dn_fermi;
            XDIAG_FILL(base_in + idx_dnc, idx_out, neg ? -cf : cf);
          }
        }
        ++idx_dnc;
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

} // namespace xdiag::kernels::tj
