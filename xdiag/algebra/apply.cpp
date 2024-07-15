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
    if (bonds.isreal()) {
      if (v.isreal() && w.isreal()) {
        arma::vec vvec = v.vector(0, false);
        arma::vec wvec = w.vector(0, false);
        apply(bonds, v.block(), vvec, w.block(), wvec);
      } else if (v.isreal() && !w.isreal()) {
        auto w2 = State(w.block(), true);
        arma::vec vvec = v.vector(0, false);
        arma::vec wvec = w2.vector(0, false);
        apply(bonds, v.block(), vvec, w.block(), wvec);
        w = w2;
      } else if (!v.isreal() && w.isreal()) {
        w.make_complex();
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(bonds, v.block(), vvec, w.block(), wvec);
      } else if (!v.isreal() && !w.isreal()) {
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(bonds, v.block(), vvec, w.block(), wvec);
      }
    } else {
      if (v.isreal()) {
        auto v2 = v;
        v2.make_complex();
        w.make_complex();
        arma::cx_vec vvec = v2.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(bonds, v.block(), vvec, w.block(), wvec);
      } else {
	w.make_complex();
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(bonds, v.block(), vvec, w.block(), wvec);
      }
    }
  } else {
    XDIAG_THROW("Applying a BondList to a state with multiple "
                "columns not yet implemented");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag
