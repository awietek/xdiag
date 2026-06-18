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

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/basis/basis_electron_symmetric.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/matrices/terms/cdagc_string.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::electron {

// Applies a product of elementary electron operators (Cdagup/Cup/Cdagdn/Cdn) to
// the electron product basis -- the spinful analogue of
// fermion::term_cdagc_string. The electron block is two independent fermion
// sectors in the "all ups then all dns" Jordan-Wigner order, so the string
// splits into an up sub-string acting on ups and a dn sub-string acting on dns,
// each handled by the shared CdagCString. The only spinful addition is the
// cross sign: every dn operator sits behind the whole up string, so the dn
// block carries (-1)^(q_dn * Nup) with q_dn the number of dn operators and Nup
// the up occupation (evaluated on the incoming ups, since in this order the dn
// operators act first).
//
// The up sub-string and the per-up cross sign are evaluated once per up state in
// the outer loop; the dn sub-string action is precomputed once per dn state.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagc_string(Coeff const &c, Monomial const &mono,
                       basis::BasisElectron<enumeration_t> const &basis_in,
                       basis::BasisElectron<enumeration_t> const &basis_out,
                       fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  // The split below evaluates the up and dn sub-strings as if the monomial were
  // in all-ups-then-all-dns order; fold in the Jordan-Wigner sign of partitioning
  // the input into that order (a no-op for the electron block's all-ups-then-dns
  // normal order, needed for the interleaved tJ creation-major normal order).
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
      XDIAG_THROW("electron Cdag/C string: unexpected operator type \"" + type +
                  "\"; only Cdagup, Cup, Cdagdn, Cdn are allowed");
    }
  }
  CdagCString<bit_t> up_str(nsites, Monomial(up_ops), "Cdagup", "Cup");
  CdagCString<bit_t> dn_str(nsites, Monomial(dn_ops), "Cdagdn", "Cdn");
  bool cross_when_odd_nup = (dn_ops.size() & 1);

  auto const &basis_up_in = basis_in.basis_up();
  auto const &basis_up_out = basis_out.basis_up();
  auto const &basis_dn_in = basis_in.basis_dn();
  auto const &basis_dn_out = basis_out.basis_dn();
  int64_t size_dn_in = basis_dn_in.size();
  int64_t size_dn_out = basis_dn_out.size();

  // Precompute the dn sub-string action for every dn state (independent of the
  // up sector). dn_out[idx] == -1 marks a dn state annihilated by the string.
  // std::vector<char> (not <bool>) avoids a data race on the packed bits.
  std::vector<int64_t> dn_out(size_dn_in);
  std::vector<char> dn_neg(size_dn_in, 0);
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_dn, end_dn, idx_dn] =
        utils::thread_range(basis_dn_in, num_thread, omp_get_num_threads());
#else
    auto [begin_dn, end_dn, idx_dn] = utils::thread_range(basis_dn_in, 0, 1);
#endif
    for (auto it = begin_dn; it != end_dn; ++it, ++idx_dn) {
      bit_t dns = *it;
      if (dn_str.non_zero(dns)) {
        std::pair<bit_t, bool> a = dn_str.action(dns);
        dn_out[idx_dn] = basis_dn_out.index(a.first);
        dn_neg[idx_dn] = a.second ? 1 : 0;
      } else {
        dn_out[idx_dn] = -1;
      }
    }
#ifdef _OPENMP
  }
#endif

  // Apply the up sub-string once per up state and fold in the dn table.
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up_in, num_thread, omp_get_num_threads());
#else
    auto [begin_up, end_up, idx_up] = utils::thread_range(basis_up_in, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      if (!up_str.non_zero(ups)) {
        continue;
      }
      std::pair<bit_t, bool> a = up_str.action(ups);
      bool up_total = a.second;
      if (cross_when_odd_nup && (bits::popcount(ups) & 1)) {
        up_total = !up_total;
      }
      int64_t base_in = idx_up * size_dn_in;
      int64_t base_out = basis_up_out.index(a.first) * size_dn_out;
      for (int64_t idx_dn = 0; idx_dn < size_dn_in; ++idx_dn) {
        if (dn_out[idx_dn] >= 0) {
          bool neg = up_total ^ (bool)dn_neg[idx_dn];
          XDIAG_FILL(base_in + idx_dn, base_out + dn_out[idx_dn],
                     neg ? -cf : cf);
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

// Symmetric overload. Same split into an up sub-string and a dn sub-string, but
// the result (ups_mid, dns_mid) is mapped into the symmetric output basis. The
// up sub-string action, the cross sign (-1)^(q_dn * Nup_in), the output up
// representative and its up-stabilizer are all evaluated once per input up
// representative (outer loop, hoisted); the inner dn loop applies the dn
// sub-string and re-symmetrises the dn within the output up rep's block. The
// total sign combines the operator signs (up, dn, cross) with the
// symmetrisation fermi sign (fermi_up XOR fermi_dn).
template <typename coeff_t, class basis_t, class fill_f,
          typename = decltype(std::declval<basis_t>().dns_for_ups_rep(0))>
void term_cdagc_string(Coeff const &c, Monomial const &mono,
                       basis_t const &basis_in, basis_t const &basis_out,
                       fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  // See the non-symmetric overload: fold in the Jordan-Wigner sign of
  // partitioning the monomial into all-ups-then-all-dns order.
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
      XDIAG_THROW("electron Cdag/C string: unexpected operator type \"" + type +
                  "\"; only Cdagup, Cup, Cdagdn, Cdn are allowed");
    }
  }
  CdagCString<bit_t> up_str(nsites, Monomial(up_ops), "Cdagup", "Cup");
  CdagCString<bit_t> dn_str(nsites, Monomial(dn_ops), "Cdagdn", "Cdn");
  bool cross_when_odd_nup = (dn_ops.size() & 1);

  auto const &basis_up_in = basis_in.basis_up();
  auto bloch = basis_out.characters().template as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up_in] =
        utils::thread_range(basis_up_in, num_thread, omp_get_num_threads());
#else
    auto [begin_up, end_up, idx_up_in] =
        utils::thread_range(basis_up_in, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up_in) {
      bit_t ups = *it_up;
      if (!up_str.non_zero(ups)) {
        continue;
      }
      // up sub-string action + cross sign, both hoisted out of the dn loop
      auto [ups_flip, up_fermi] = up_str.action(ups);
      if (cross_when_odd_nup && (bits::popcount(ups) & 1)) {
        up_fermi = !up_fermi;
      }

      auto [raw, s0, nrm] = basis_out.basis_up().representative_data(ups_flip);
      (void)nrm;
      int64_t idx_up_out = raw - 1;
      int64_t off_in = basis_in.ups_offset(idx_up_in);
      int64_t off_out = basis_out.ups_offset(idx_up_out);
      auto dnss_in = basis_in.dns_for_ups_rep(idx_up_in);
      auto norms_in = basis_in.norms_for_ups_rep(idx_up_in);

      if (basis_out.stab_size(idx_up_out) == 1) {
        coeff_t prefac = cf * bloch(s0);
        bool fermi_up = basis_out.fermi_bool_ups(s0, ups_flip);
        auto dnss_out = basis_out.dns_for_ups_rep(idx_up_out);
        int64_t dn_idx = 0;
        for (bit_t dns : dnss_in) {
          if (dn_str.non_zero(dns)) {
            auto [dns_flip, dn_fermi] = dn_str.action(dns);
            auto [idx_dns_out, fermi_dn] =
                basis_out.index_dns_fermi(dns_flip, s0, idx_up_out, dnss_out);
            if (idx_dns_out >= 0) { // -1: tJ double occupancy (never for electron)
              coeff_t val = prefac / norms_in[dn_idx];
              bool neg = up_fermi ^ dn_fermi ^ fermi_up ^ fermi_dn;
              XDIAG_FILL(off_in + dn_idx, off_out + idx_dns_out,
                         neg ? -val : val);
            }
          }
          ++dn_idx;
        }
      } else {
        std::vector<int64_t> syms = basis_out.syms_ups(ups_flip);
        auto dnss_out = basis_out.dns_for_ups_rep(idx_up_out);
        auto norms_out = basis_out.norms_for_ups_rep(idx_up_out);
        int64_t dn_idx = 0;
        for (bit_t dns : dnss_in) {
          if (dn_str.non_zero(dns)) {
            auto [dns_flip, dn_fermi] = dn_str.action(dns);
            auto [idx_dns_out, fermi_dn, sym] =
                basis_out.index_dns_fermi_sym(dns_flip, syms, dnss_out);
            if (idx_dns_out >= 0) {
              bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
              coeff_t val = cf * bloch(sym) * norms_out[idx_dns_out] /
                            norms_in[dn_idx];
              bool neg = up_fermi ^ dn_fermi ^ fermi_up ^ fermi_dn;
              XDIAG_FILL(off_in + dn_idx, off_out + idx_dns_out,
                         neg ? -val : val);
            }
          }
          ++dn_idx;
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
