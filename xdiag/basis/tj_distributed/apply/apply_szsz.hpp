#pragma once

#include <xdiag/basis/tj_distributed/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj_distributed {

template <typename bit_t, typename coeff_t, class basis_t>
void apply_szsz(Coupling const &cpl, Op const &op, basis_t const &basis,
                const coeff_t *vec_in, coeff_t *vec_out) {
  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  std::string type = op.type();
  bit_t s1_mask = (bit_t)1 << s1;
  bit_t s2_mask = (bit_t)1 << s2;

  // Set values for same/diff (tJ block definition)
  coeff_t val_same, val_diff;
  if (type == "SzSz") {
    val_same = J / 4.;
    val_diff = -J / 4.;
  } else { // (type == "tJSzSz")
    val_same = 0.;
    val_diff = -J / 2.;
  }

  auto term_action = [&](bit_t up, bit_t dn) -> coeff_t {
    bool b1_up = up & s1_mask;
    bool b2_up = up & s2_mask;
    bool b1_dn = dn & s1_mask;
    bool b2_dn = dn & s2_mask;

    if ((b1_up && b2_up) || (b1_dn && b2_dn)) {
      return val_same;
    } else if ((b1_up && b2_dn) || (b1_dn && b2_up)) {
      return val_diff;
    } else {
      return 0.;
    }
  };

  tj_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action, vec_in,
                                                    vec_out);
}

} // namespace xdiag::basis::tj_distributed
