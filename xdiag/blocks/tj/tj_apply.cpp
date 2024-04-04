#include "tj_apply.h"

#include <xdiag/algebra/fill.h>
#include <xdiag/blocks/tj/compile.h>
#include <xdiag/blocks/tj/dispatch.h>

namespace xdiag {
template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out) try {
  vec_out.zeros();

  BondList bondsc = tj::compile(bonds, 1e-12);
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  tj::dispatch<coeff_t>(bondsc, block_in, block_out, fill);
} catch (...) {
  XDiagRethrow("Cannot apply bonds on vector for \"tJ\" block");
}

template void apply<double>(BondList const &, tJ const &,
                            arma::Col<double> const &, tJ const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, tJ const &,
                             arma::Col<complex> const &, tJ const &,
                             arma::Col<complex> &);

} // namespace xdiag
