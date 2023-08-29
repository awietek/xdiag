#include "tj_apply.h"

#include <hydra/algebra/generic_operator.h>
#include <hydra/blocks/tj/compile.h>

namespace hydra {
template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out) {
  generic_apply(bonds, block_in, vec_in, block_out, vec_out, tj::compile);
}

template void apply<double>(BondList const &, tJ const &,
                            arma::Col<double> const &, tJ const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, tJ const &,
                             arma::Col<complex> const &, tJ const &,
                             arma::Col<complex> &);

} // namespace hydra
