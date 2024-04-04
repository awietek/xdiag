#pragma once

#include <xdiag/common.h>
#include <xdiag/utils/print_macro.h>

#include <xdiag/blocks/electron/compile.h>

#include <xdiag/blocks/electron/terms/apply_exchange.h>
#include <xdiag/blocks/electron/terms/apply_hopping.h>
#include <xdiag/blocks/electron/terms/apply_ising.h>
#include <xdiag/blocks/electron/terms/apply_number.h>
#include <xdiag/blocks/electron/terms/apply_raise_lower.h>
#include <xdiag/blocks/electron/terms/apply_u.h>

namespace xdiag::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(BondList const &bonds, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) {

  BondList bonds_compiled = electron::compile(bonds);

  // Separate bond acting on ups, dns, diagonally or fully mixing
  BondList bonds_ups;
  BondList bonds_dns;
  BondList bonds_diag;
  BondList bonds_mixed;
  for (auto bond : bonds_compiled) {
    std::string type = bond.type();
    if ((type == "HOPUP") || (type == "CDAGUP") || (type == "CUP")) {
      bonds_ups << bond;
    } else if ((type == "HOPDN") || (type == "CDAGDN") || (type == "CDN")) {
      bonds_dns << bond;
    } else if ((type == "ISING") || (type == "NUMBERUP") ||
               (type == "NUMBERDN")) {
      bonds_diag << bond;
    } else if (type == "EXCHANGE") {
      bonds_mixed << bond;
    } else {
      XDiagPrint(bond);
      Log.err("Error: Unknown bond of type {}", type);
    }
  }

  // Diagonal terms
  for (auto bond : bonds_diag) {
    std::string type = bond.type();
    if (type == "ISING") {
      electron::apply_ising<bit_t, coeff_t, symmetric>(bond, basis_in, fill);
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      electron::apply_number<bit_t, coeff_t, symmetric>(bond, basis_in, fill);
    }
  }

  if (bonds.coupling_defined("U")) {
    coeff_t U = bonds.coupling<coeff_t>("U");
    electron::apply_u<bit_t, coeff_t, symmetric>(U, basis_in, fill);
  }

  // terms on both ups and dns
  for (auto bond : bonds_mixed) {
    std::string type = bond.type();
    if (type == "EXCHANGE") {
      electron::apply_exchange<bit_t, coeff_t, symmetric>(bond, basis_in, fill);
    }
  }

  // terms acting only on ups
  for (auto bond : bonds_ups) {
    std::string type = bond.type();
    if (type == "HOPUP") {
      electron::apply_hopping<bit_t, coeff_t, symmetric>(bond, basis_in, fill);
    } else if ((type == "CDAGUP") || (type == "CUP")) {
      electron::apply_raise_lower<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                             basis_out, fill);
    }
  }

  // terms acting only on dns
  for (auto bond : bonds_dns) {
    std::string type = bond.type();
    if (type == "HOPDN") {
      electron::apply_hopping<bit_t, coeff_t, symmetric>(bond, basis_in, fill);
    } else if ((type == "CDAGDN") || (type == "CDN")) {
      electron::apply_raise_lower<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                             basis_out, fill);
    }
  }
}

} // namespace xdiag::electron
