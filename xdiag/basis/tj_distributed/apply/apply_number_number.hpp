#pragma once

#include <xdiag/basis/tj_distributed/apply/generic_term_diag.hpp>
#include <xdiag/bits/gbit.hpp>

namespace xdiag::basis::tj_distributed {

template <typename coeff_t, class basis_t>
void apply_ntot_ntot(Coupling const &cpl, Op const &op, basis_t const &basis,
                     const coeff_t *vec_in, coeff_t *vec_out) try {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  auto apply = [&](bit_t ups, bit_t dns) {
    int n1 = gbit(ups, s1) + gbit(dns, s1);
    int n2 = gbit(ups, s2) + gbit(dns, s2);
    return mu * (coeff_t)(n1 * n2);
  };

  tj_distributed::generic_term_diag<coeff_t>(basis, apply, vec_in, vec_out);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj_distributed
