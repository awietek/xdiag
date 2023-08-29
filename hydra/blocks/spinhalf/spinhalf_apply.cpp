#include "spinhalf_apply.h"

#include <hydra/algebra/generic_operator.h>
#include <hydra/blocks/spinhalf/compile.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, Spinhalf const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
           arma::Col<coeff_t> &vec_out) {
  generic_apply(bonds, block_in, vec_in, block_out, vec_out, spinhalf::compile);
}

template void apply<double>(BondList const &, Spinhalf const &,
                            arma::Col<double> const &, Spinhalf const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, Spinhalf const &,
                             arma::Col<complex> const &, Spinhalf const &,
                             arma::Col<complex> &);

} // namespace hydra
