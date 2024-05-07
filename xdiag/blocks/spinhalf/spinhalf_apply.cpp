#include "spinhalf_apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/spinhalf/compile.hpp>
#include <xdiag/blocks/spinhalf/dispatch.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(BondList const &bonds, Spinhalf const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
           arma::Col<coeff_t> &vec_out) try {
  vec_out.zeros();

  BondList bondsc = spinhalf::compile(bonds, 1e-12);
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  spinhalf::dispatch<coeff_t>(bondsc, block_in, block_out, fill);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(BondList const &, Spinhalf const &,
                            arma::Col<double> const &, Spinhalf const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, Spinhalf const &,
                             arma::Col<complex> const &, Spinhalf const &,
                             arma::Col<complex> &);

} // namespace xdiag
