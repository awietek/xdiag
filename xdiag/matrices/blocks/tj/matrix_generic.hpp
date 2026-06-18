// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <type_traits>
#include <utility>

#include <xdiag/algebra/algebras/tj_implementation_algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

#include <xdiag/matrices/blocks/tj/terms/term_cdagc_string.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_exchange.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_hopdn.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_hopup.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_number.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_raise_lower.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_szsz.hpp>
#include <xdiag/matrices/terms/term_identity.hpp>

// The permutation-symmetric tJ basis (BasistJSymmetric) reuses the electron
// symmetric matrix kernels (templated on the basis type); see matrix_generic
// below for the dispatch.
#include <xdiag/matrices/blocks/electron/terms/term_cdagc_string.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_diag.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_hop.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_number.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_number_number.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_raise_lower.hpp>

namespace xdiag::matrices::tj {

// Detects the coupled-symmetric tJ basis (it exposes dns_for_ups_rep, the
// compressed non-symmetric BasistJ does not).
template <typename T, typename = void>
struct is_symmetric_basis : std::false_type {};
template <typename T>
struct is_symmetric_basis<
    T, std::void_t<decltype(std::declval<T>().dns_for_ups_rep(0))>>
    : std::true_type {};

// Symmetric tJ matrix kernel: same tJ algebra (projected operators, tJ
// normal-ordering, no double occupancy) as the non-symmetric block, but the
// per-term work is delegated to the shared electron SYMMETRIC kernels, which
// drop any double-occupied output via the basis' index_dns returning -1. Only
// two pieces are tJ-specific: the tJSzSz diagonal value (electron has SzSz only)
// and Exchange (routed through the electron Cdag/C string kernel).
template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic_symmetric(OpSum const &ops, basis_t const &basis_in,
                              basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  operators::check_valid(ops);
  // No Exchange kernel for the symmetric basis: let it expand to (normal-ordered)
  // Cdag/C strings handled by the electron string kernel.
  auto algebra = algebra::tj_implementation_algebra(basis_in.nsites(),
                                                    /*exchange_as_kernel=*/false);
  auto ops_compiled = normal_order(ops.plain(), algebra);

  for (auto const &[c, monomial] : ops_compiled) {
    if (monomial.size() != 1) {
      electron::term_cdagc_string<coeff_t>(c, monomial, basis_in, basis_out,
                                           fill);
      continue;
    }
    Op op = monomial[0];
    std::string type = op.type();

    if (type == "Id") {
      matrices::term_identity<coeff_t>(c, basis_in, fill);
    } else if (type == "Hopup") {
      electron::term_hopup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Hopdn") {
      electron::term_hopdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cup") {
      electron::term_cup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagup") {
      electron::term_cdagup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdn") {
      electron::term_cdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagdn") {
      electron::term_cdagdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Nup") {
      electron::term_nup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Ndn") {
      electron::term_ndn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NupNup") {
      electron::term_nupnup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NdnNdn") {
      electron::term_ndnndn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NupNdn") {
      electron::term_nupndn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NdnNup") {
      electron::term_ndnnup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NtotNtot") {
      electron::term_ntotntot<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "SzSz") || (type == "tJSzSz")) {
      // diagonal spin-spin coupling: SzSz = +/-J/4, tJSzSz = 0 / -J/2 (see the
      // non-symmetric tj::term_szsz for the convention).
      coeff_t J = c.scalar().template as<coeff_t>();
      int64_t i = op[0], j = op[1];
      coeff_t val_same =
          (type == "tJSzSz") ? coeff_t(0.0) : J / coeff_t(4.0);
      coeff_t val_diff =
          (type == "tJSzSz") ? -J / coeff_t(2.0) : -J / coeff_t(4.0);
      electron::term_diag<coeff_t>(
          basis_in, basis_out,
          [=](bit_t ups, bit_t dns) -> coeff_t {
            bool up_i = bits::get(ups, i), up_j = bits::get(ups, j);
            bool dn_i = bits::get(dns, i), dn_j = bits::get(dns, j);
            if ((up_i && up_j) || (dn_i && dn_j)) {
              return val_same;
            } else if ((up_i && dn_j) || (dn_i && up_j)) {
              return val_diff;
            }
            return coeff_t(0.0);
          },
          fill);
    } else {
      XDIAG_THROW(
          fmt::format("Unknown Op type for symmetric tJ basis \"{}\"", type));
    }
  }
}
XDIAG_CATCH

template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic(OpSum const &ops, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) try {

  if constexpr (is_symmetric_basis<basis_t>::value) {
    matrix_generic_symmetric<coeff_t>(ops, basis_in, basis_out, fill);
  } else {

  operators::check_valid(ops);
  auto algebra = algebra::tj_implementation_algebra(basis_in.nsites());
  auto ops_compiled = normal_order(ops.plain(), algebra);

  for (auto const &[c, monomial] : ops_compiled) {
    // Products of elementary operators are handled by the general Cdag/C string
    // kernel; size-1 named operators dispatch to their dedicated fast kernels.
    if (monomial.size() != 1) {
      term_cdagc_string<coeff_t>(c, monomial, basis_in, basis_out, fill);
      continue;
    }
    Op op = monomial[0];
    std::string type = op.type();

    if (type == "Id") {
      matrices::term_identity<coeff_t>(c, basis_in, fill);
    } else if (type == "Hopup") {
      term_hopup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Hopdn") {
      term_hopdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cup") {
      term_cup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagup") {
      term_cdagup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdn") {
      term_cdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagdn") {
      term_cdagdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "Nup") || (type == "Ndn")) {
      term_number<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "NupNup") || (type == "NdnNdn") || (type == "NupNdn") ||
               (type == "NdnNup") || (type == "NtotNtot")) {
      term_number_number<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "SzSz") || (type == "tJSzSz")) {
      term_szsz<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Exchange") {
      term_exchange<coeff_t>(c, op, basis_in, basis_out, fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type for tJ basis \"{}\"", type));
    }
  }
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices::tj
