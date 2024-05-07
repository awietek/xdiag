#include "apply.hpp"

#include <variant>
#include <xdiag/common.hpp>

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
} catch (Error const &error) {
  XDIAG_RETHROW(error);
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
      XDIAG_THROW(
          "Apply operator only works if both states are complex or both states "
          "are real and the operator is real. Consider Making the vectors "
          "complex first by using method .make_complex().");
    }
  } else {
    XDIAG_THROW("Applying a BondList to a state with multiple "
                "columns not yet implemented");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag
