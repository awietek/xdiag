#include "electron_apply.h"

#include <hydra/blocks/electron/terms/electron_exchange.h>
#include <hydra/blocks/electron/terms/electron_hopping.h>
#include <hydra/blocks/electron/terms/electron_ising.h>
#include <hydra/blocks/electron/terms/electron_u.h>

#include <hydra/blocks/electron/terms/electron_symmetric_exchange.h>
#include <hydra/blocks/electron/terms/electron_symmetric_hopping.h>
#include <hydra/blocks/electron/terms/electron_symmetric_ising.h>
#include <hydra/blocks/electron/terms/electron_symmetric_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, arma::Col<coeff_t> const &vec_in,
           Electron<bit_t> const &block_out, arma::Col<coeff_t> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == (idx_t)vec_in.size());
  assert(block_out.size() == (idx_t)vec_out.size());

  utils::check_operator_works_with<coeff_t>(bonds, couplings, "tj_apply");

  vec_out.zeros();
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  if (block_in.symmetric()) {

    if (block_in.charge_conserved() && block_in.sz_conserved()) {
      auto const &indexing_in = block_in.indexing_sym_np();
      terms::electron_symmetric_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_symmetric_hopping<bit_t, coeff_t>(bonds, couplings,
                                                        indexing_in, fill);
      terms::electron_symmetric_ising<bit_t, coeff_t>(bonds, couplings,
                                                      indexing_in, fill);
      terms::electron_symmetric_exchange<bit_t, coeff_t>(bonds, couplings,
                                                         indexing_in, fill);
    } else {
      auto const &indexing_in = block_in.indexing_sym_no_np();
      terms::electron_symmetric_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_symmetric_hopping<bit_t, coeff_t>(bonds, couplings,
                                                        indexing_in, fill);
      terms::electron_symmetric_ising<bit_t, coeff_t>(bonds, couplings,
                                                      indexing_in, fill);
      terms::electron_symmetric_exchange<bit_t, coeff_t>(bonds, couplings,
                                                         indexing_in, fill);
    }

  } else {

    if (block_in.charge_conserved() && block_in.sz_conserved()) {
      auto const &indexing_in = block_in.indexing_np();
      terms::electron_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                              fill);
      terms::electron_ising<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                            fill);
      terms::electron_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                               fill);
    } else {
      auto const &indexing_in = block_in.indexing_no_np();
      terms::electron_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                              fill);
      terms::electron_ising<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                            fill);
      terms::electron_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                               fill);
    }
  }
}

template void Apply<uint16_t, double>(BondList const &, Couplings const &,
                                      Electron<uint16_t> const &,
                                      arma::Col<double> const &,
                                      Electron<uint16_t> const &,
                                      arma::Col<double> &);
template void Apply<uint32_t, double>(BondList const &, Couplings const &,
                                      Electron<uint32_t> const &,
                                      arma::Col<double> const &,
                                      Electron<uint32_t> const &,
                                      arma::Col<double> &);
template void Apply<uint64_t, double>(BondList const &, Couplings const &,
                                      Electron<uint64_t> const &,
                                      arma::Col<double> const &,
                                      Electron<uint64_t> const &,
                                      arma::Col<double> &);

template void Apply<uint16_t, complex>(BondList const &, Couplings const &,
                                       Electron<uint16_t> const &,
                                       arma::Col<complex> const &,
                                       Electron<uint16_t> const &,
                                       arma::Col<complex> &);
template void Apply<uint32_t, complex>(BondList const &, Couplings const &,
                                       Electron<uint32_t> const &,
                                       arma::Col<complex> const &,
                                       Electron<uint32_t> const &,
                                       arma::Col<complex> &);
template void Apply<uint64_t, complex>(BondList const &, Couplings const &,
                                       Electron<uint64_t> const &,
                                       arma::Col<complex> const &,
                                       Electron<uint64_t> const &,
                                       arma::Col<complex> &);

} // namespace hydra
