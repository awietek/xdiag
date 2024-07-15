#include "tj_apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/tj/compile.hpp>
#include <xdiag/blocks/tj/dispatch.hpp>

namespace xdiag {
template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision) try {
  vec_out.zeros();

  BondList bondsc = tj::compile(bonds, block_in.n_sites(), zero_precision);
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  tj::dispatch<coeff_t>(bondsc, block_in, block_out, fill, zero_precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(BondList const &, tJ const &,
                            arma::Col<double> const &, tJ const &,
                            arma::Col<double> &, double);

template void apply<complex>(BondList const &, tJ const &,
                             arma::Col<complex> const &, tJ const &,
                             arma::Col<complex> &, double);

} // namespace xdiag
