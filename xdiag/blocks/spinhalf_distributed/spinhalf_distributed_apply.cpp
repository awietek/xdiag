#include "tj_apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/blocks/tj/compile.hpp>
#include <xdiag/blocks/tj/dispatch.hpp>

namespace xdiag {
template <typename coeff_t>
void apply(OpSum const &ops, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision) try {
  vec_out.zeros();

  OpSum opsc = spinhalf::compile(ops, block_in, zero_precision);
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  tj::dispatch<coeff_t>(opsc, block_in, block_out, fill);
} catch (...) {
  XDiagRethrow("Cannot apply ops on vector for \"tJ\" block");
}

template void apply<double>(OpSum const &, tJ const &,
                            arma::Col<double> const &, tJ const &,
                            arma::Col<double> &, double);

template void apply<complex>(OpSum const &, tJ const &,
                             arma::Col<complex> const &, tJ const &,
                             arma::Col<complex> &, double);

} // namespace xdiag
