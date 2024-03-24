#include "tj_apply.h"

#include <hydra/algebra/fill.h>
#include <hydra/blocks/tj/compile.h>
#include <hydra/blocks/tj/dispatch.h>

namespace hydra {
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
  HydraRethrow("Cannot apply bonds on vector for \"tJ\" block");
}

template void apply<double>(BondList const &, tJ const &,
                            arma::Col<double> const &, tJ const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, tJ const &,
                             arma::Col<complex> const &, tJ const &,
                             arma::Col<complex> &);

} // namespace hydra
