// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_diag.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels::tj {

// Single-site number operators Nup{i} / Ndn{i} (diagonal). On the tJ basis a
// site is empty / up / dn, so n^up and n^dn are 0/1 and n^up n^dn == 0.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_number(Coeff const &c, Op const &op,
                 basis::BasistJ<enumeration_t> const &basis_in,
                 basis::BasistJ<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t i = op[0];
  int64_t nsites = basis_in.nsites();
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);
  bool is_up = (op.type() == "Nup");

  term_diag<coeff_t>(
      basis_in, basis_out,
      [=](bit_t ups) {
        bool up_i = bits::get(ups, i);
        int64_t ri = tj_dn_rank<bit_t>(ups, i, nsites, sitesmask);
        return [=](bit_t dnc) -> coeff_t {
          bool occ = is_up ? up_i : ((ri >= 0) && bits::get(dnc, ri));
          return occ ? cf : coeff_t(0.0);
        };
      },
      fill);
}
XDIAG_CATCH

// Two-site number products on a bond {i,j} (diagonal): NtotNtot, NupNup, NdnNdn,
// NupNdn, NdnNup. Ntot_s = n^up_s + n^dn_s (= occupancy, 0/1 on the tJ basis).
template <typename coeff_t, class enumeration_t, class fill_f>
void term_number_number(Coeff const &c, Op const &op,
                        basis::BasistJ<enumeration_t> const &basis_in,
                        basis::BasistJ<enumeration_t> const &basis_out,
                        fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t i = op[0];
  int64_t j = op[1];
  int64_t nsites = basis_in.nsites();
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);

  // Resolve the type ONCE (not in the hot loop) to which sectors enter each
  // factor: n_i = (want_up_i ? n^up_i) | (want_dn_i ? n^dn_i), likewise n_j.
  std::string type = op.type();
  bool wui, wdi, wuj, wdj;
  if (type == "NupNup") {
    wui = true, wdi = false, wuj = true, wdj = false;
  } else if (type == "NdnNdn") {
    wui = false, wdi = true, wuj = false, wdj = true;
  } else if (type == "NupNdn") {
    wui = true, wdi = false, wuj = false, wdj = true;
  } else if (type == "NdnNup") {
    wui = false, wdi = true, wuj = true, wdj = false;
  } else { // NtotNtot
    wui = true, wdi = true, wuj = true, wdj = true;
  }

  term_diag<coeff_t>(
      basis_in, basis_out,
      [=](bit_t ups) {
        // up contributions are fixed per ups; only the dn parts vary with dnc
        bool ni_up = wui && bits::get(ups, i);
        bool nj_up = wuj && bits::get(ups, j);
        int64_t ri = tj_dn_rank<bit_t>(ups, i, nsites, sitesmask);
        int64_t rj = tj_dn_rank<bit_t>(ups, j, nsites, sitesmask);
        return [=](bit_t dnc) -> coeff_t {
          bool ni = ni_up || (wdi && (ri >= 0) && bits::get(dnc, ri));
          bool nj = nj_up || (wdj && (rj >= 0) && bits::get(dnc, rj));
          return (ni && nj) ? cf : coeff_t(0.0);
        };
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::kernels::tj
