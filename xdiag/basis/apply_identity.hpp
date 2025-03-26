#pragma once
#include <xdiag/common.hpp>
#include <xdiag/operators/coupling.hpp>

namespace xdiag::basis {

template <typename coeff_t, class basis_t, class fill_f>
void apply_identity(Coupling const &cpl, basis_t const &basis,
                    fill_f fill) try {
  coeff_t s = cpl.scalar().as<coeff_t>();
  for (int64_t i = 0; i < basis.size(); ++i) {
    fill(i, i, s);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag::basis
