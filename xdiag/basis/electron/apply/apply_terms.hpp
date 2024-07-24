#pragma once

#include <xdiag/basis/electron/apply/apply_exchange.hpp>
#include <xdiag/basis/electron/apply/apply_hopping.hpp>
#include <xdiag/basis/electron/apply/apply_ising.hpp>
#include <xdiag/basis/electron/apply/apply_number.hpp>
#include <xdiag/basis/electron/apply/apply_raise_lower.hpp>
#include <xdiag/basis/electron/apply/apply_u.hpp>

#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) try {
  // Separate op acting on ups, dns, diagonally or fully mixing
  OpSum ops_ups;
  OpSum ops_dns;
  OpSum ops_diag;
  OpSum ops_mixed;
  for (auto op : ops) {
    std::string type = op.type();
    if ((type == "HOPUP") || (type == "CDAGUP") || (type == "CUP")) {
      ops_ups += op;
    } else if ((type == "HOPDN") || (type == "CDAGDN") || (type == "CDN")) {
      ops_dns += op;
    } else if ((type == "ISING") || (type == "NUMBERUP") ||
               (type == "NUMBERDN")) {
      ops_diag += op;
    } else if (type == "EXCHANGE") {
      ops_mixed += op;
    } else {
      XDIAG_THROW(fmt::format("Error: Unknown Op of type \"{}\"", type));
    }
  }

  // Diagonal terms
  for (auto op : ops_diag) {
    std::string type = op.type();
    if (type == "ISING") {
      electron::apply_ising<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      electron::apply_number<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    }
  }

  if (ops.defined("U")) {
    Coupling cpl = ops["U"];

    if (cpl.is<double>()) {
      double U = cpl.as<double>();
      electron::apply_u<bit_t, coeff_t, symmetric>(U, basis_in, fill);
    } else {
      XDIAG_THROW("Coupling U must either be a real number");
    }
  }

  // terms on both ups and dns
  for (auto op : ops_mixed) {
    std::string type = op.type();
    if (type == "EXCHANGE") {
      electron::apply_exchange<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    }
  }

  // terms acting only on ups
  for (auto op : ops_ups) {
    std::string type = op.type();
    if (type == "HOPUP") {
      electron::apply_hopping<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if ((type == "CDAGUP") || (type == "CUP")) {
      electron::apply_raise_lower<bit_t, coeff_t, symmetric>(op, basis_in,
                                                             basis_out, fill);
    }
  }

  // terms acting only on dns
  for (auto op : ops_dns) {
    std::string type = op.type();
    if (type == "HOPDN") {
      electron::apply_hopping<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if ((type == "CDAGDN") || (type == "CDN")) {
      electron::apply_raise_lower<bit_t, coeff_t, symmetric>(op, basis_in,
                                                             basis_out, fill);
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
