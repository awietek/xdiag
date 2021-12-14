#include "electron_symmetric_apply.h"

#include <hydra/blocks/electron_symmetric/electron_symmetric.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_exchange.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_hopping.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_ising.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           ElectronSymmetric<bit_t> const &block_out,
           lila::Vector<coeff_t> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_operator_works_with<coeff_t>(bonds, couplings, block_in.irrep(),
                                            block_out.irrep(),
                                            "electron_symmetric_apply");
  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  using namespace terms;
  electron_symmetric_U<bit_t, coeff_t>(couplings, indexing_in, fill);
  electron_symmetric_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                             fill);
  electron_symmetric_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  electron_symmetric_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                              fill);
}

template void Apply<uint16_t, double>(BondList const &, Couplings const &,
                                      ElectronSymmetric<uint16_t> const &,
                                      lila::Vector<double> const &,
                                      ElectronSymmetric<uint16_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint32_t, double>(BondList const &, Couplings const &,
                                      ElectronSymmetric<uint32_t> const &,
                                      lila::Vector<double> const &,
                                      ElectronSymmetric<uint32_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint64_t, double>(BondList const &, Couplings const &,
                                      ElectronSymmetric<uint64_t> const &,
                                      lila::Vector<double> const &,
                                      ElectronSymmetric<uint64_t> const &,
                                      lila::Vector<double> &);

template void Apply<uint16_t, complex>(BondList const &, Couplings const &,
                                       ElectronSymmetric<uint16_t> const &,
                                       lila::Vector<complex> const &,
                                       ElectronSymmetric<uint16_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint32_t, complex>(BondList const &, Couplings const &,
                                       ElectronSymmetric<uint32_t> const &,
                                       lila::Vector<complex> const &,
                                       ElectronSymmetric<uint32_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint64_t, complex>(BondList const &, Couplings const &,
                                       ElectronSymmetric<uint64_t> const &,
                                       lila::Vector<complex> const &,
                                       ElectronSymmetric<uint64_t> const &,
                                       lila::Vector<complex> &);

} // namespace hydra
