#include "tj_apply.h"

#include <hydra/blocks/tj/terms/tj_exchange.h>
#include <hydra/blocks/tj/terms/tj_hopping.h>
#include <hydra/blocks/tj/terms/tj_ising.h>

#include <hydra/blocks/tj/terms/tj_symmetric_exchange.h>
#include <hydra/blocks/tj/terms/tj_symmetric_hopping.h>
#include <hydra/blocks/tj/terms/tj_symmetric_ising.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJ<bit_t> const &block_in, lila::Vector<coeff_t> const &vec_in,
           tJ<bit_t> const &block_out, lila::Vector<coeff_t> &vec_out) {
  using namespace terms;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_operator_works_with<coeff_t>(bonds, couplings, "tj_apply");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  if (block_in.symmetric()) {

    if (block_in.charge_conserved() && block_in.sz_conserved()) {
      auto const &indexing_in = block_in.indexing_sym_np();
      tj_symmetric_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
      tj_symmetric_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
      tj_symmetric_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                            fill);
    }

  } else {

    if (block_in.charge_conserved() && block_in.sz_conserved()) {
      auto const &indexing_in = block_in.indexing_np();
      tj_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
      tj_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
      tj_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
    }
  }
}

template void Apply<uint16_t, double>(BondList const &, Couplings const &,
                                      tJ<uint16_t> const &,
                                      lila::Vector<double> const &,
                                      tJ<uint16_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint32_t, double>(BondList const &, Couplings const &,
                                      tJ<uint32_t> const &,
                                      lila::Vector<double> const &,
                                      tJ<uint32_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint64_t, double>(BondList const &, Couplings const &,
                                      tJ<uint64_t> const &,
                                      lila::Vector<double> const &,
                                      tJ<uint64_t> const &,
                                      lila::Vector<double> &);

template void Apply<uint16_t, complex>(BondList const &, Couplings const &,
                                       tJ<uint16_t> const &,
                                       lila::Vector<complex> const &,
                                       tJ<uint16_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint32_t, complex>(BondList const &, Couplings const &,
                                       tJ<uint32_t> const &,
                                       lila::Vector<complex> const &,
                                       tJ<uint32_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint64_t, complex>(BondList const &, Couplings const &,
                                       tJ<uint64_t> const &,
                                       lila::Vector<complex> const &,
                                       tJ<uint64_t> const &,
                                       lila::Vector<complex> &);

} // namespace hydra
