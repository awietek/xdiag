#include "tj_distributed_apply.h"

#include <hydra/algebra/fill.h>
#include <hydra/blocks/tj/compile.h>
#include <hydra/blocks/tj_distributed/dispatch.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out) try {
  vec_out.zeros();

  BondList bondsc = tj::compile(bonds, 1e-12);
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  tj_distributed::dispatch<coeff_t>(bondsc, block_in, block_out, fill);
} catch (...) {
  HydraRethrow("Cannot apply bonds on vector for \"tJDistributed\" block");
}

template void apply<double>(BondList const &, tJDistributed const &,
                            arma::Col<double> const &, tJDistributed const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, tJDistributed const &,
                             arma::Col<complex> const &, tJDistributed const &,
                             arma::Col<complex> &);

} // namespace hydra
