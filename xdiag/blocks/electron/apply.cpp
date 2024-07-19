#include "apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/electron/compile.hpp>
#include <xdiag/blocks/electron/dispatch.hpp>

namespace xdiag {

template <typename coeff_t>
void apply(OpSum const &ops, Electron const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision) try {
  vec_out.zeros();
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  electron::dispatch<coeff_t>(ops, block_in, block_out, fill, zero_precision);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, Electron const &,
                            arma::Col<double> const &, Electron const &,
                            arma::Col<double> &, double);

template void apply<complex>(OpSum const &, Electron const &,
                             arma::Col<complex> const &, Electron const &,
                             arma::Col<complex> &, double);

} // namespace xdiag
