#include "apply.hpp"

#include <xdiag/common.hpp>
#include <variant>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, block_variant_t const &block_in,
           arma::Col<coeff_t> const &vec_in, block_variant_t const &block_out,
           arma::Col<coeff_t> &vec_out) try {
  std::visit(
      [&](auto &&block_in, auto &&block_out) {
        apply(bonds, block_in, vec_in, block_out, vec_out);
      },
      block_in, block_out);
} catch (...) {
  XDiagRethrow("Error dispatching apply");
}

template void apply(BondList const &, block_variant_t const &,
                    arma::vec const &, block_variant_t const &, arma::vec &);
template void apply(BondList const &, block_variant_t const &,
                    arma::cx_vec const &, block_variant_t const &,
                    arma::cx_vec &);

void apply(BondList const &bonds, State const &v, State &w) try {
  if ((v.n_cols() == 1) && (w.n_cols() == 1)) {
    if (bonds.isreal() && v.isreal() && w.isreal()) {
      arma::vec vvec = v.vector(0, false);
      arma::vec wvec = w.vector(0, false);
      apply(bonds, v.block(), vvec, w.block(), wvec);
    } else if (v.iscomplex() && w.iscomplex()) {
      arma::cx_vec vvec = v.vectorC(0, false);
      arma::cx_vec wvec = w.vectorC(0, false);
      apply(bonds, v.block(), vvec, w.block(), wvec);
    } else {
      XDiagThrow(
          std::logic_error,
          "Apply operator only works if both states are complex or both states "
          "are real and the operator is real. Consider Making the vectors "
          "complex first by using method .make_complex().");
    }
  } else {
    XDiagThrow(std::runtime_error,
               "Applying a BondList to a state with multiple "
               "columns not yet implemented");
  }
} catch (...) {
  XDiagRethrow("Cannot apply BondList to State");
}

} // namespace xdiag
