#pragma once

#include <xdiag/basis/tj/apply/apply_exchange.hpp>
#include <xdiag/basis/tj/apply/apply_hopping.hpp>
#include <xdiag/basis/tj/apply/apply_number.hpp>
#include <xdiag/basis/tj/apply/apply_number_number.hpp>
#include <xdiag/basis/tj/apply/apply_raise_lower.hpp>
#include <xdiag/basis/tj/apply/apply_szsz.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class basis_t,
          class fill_f>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 basis_t const &basis_out, fill_f fill) try {

  for (auto const &[cpl, op] : ops) {
    std::string type = op.type();
    if ((type == "SzSz") || (type == "tJSzSz")) {
      tj::apply_szsz<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "Nup") || (type == "Ndn")) {
      tj::apply_number<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if (type == "NtotNtot") {
      tj::apply_number_number<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                         fill);
    } else if (type == "Exchange") {
      tj::apply_exchange<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "Hopup") || (type == "Hopdn")) {
      tj::apply_hopping<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "Cdagup") || (type == "Cup") || (type == "Cdagdn") ||
               (type == "Cdn")) {
      tj::apply_raise_lower<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                       basis_out, fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type for tJ block: \"{}\"", type));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj
