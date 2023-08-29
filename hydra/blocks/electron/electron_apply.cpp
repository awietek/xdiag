#include "electron_apply.h"

#include <hydra/algebra/generic_operator.h>
#include <hydra/blocks/electron/compile.h>
#include <hydra/operators/compiler.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, Electron const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron const &block_out,
           arma::Col<coeff_t> &vec_out) {
  generic_apply(bonds, block_in, vec_in, block_out, vec_out, electron::compile);
}

template void apply<double>(BondList const &, Electron const &,
                            arma::Col<double> const &, Electron const &,
                            arma::Col<double> &);

template void apply<complex>(BondList const &, Electron const &,
                             arma::Col<complex> const &, Electron const &,
                             arma::Col<complex> &);

} // namespace hydra
