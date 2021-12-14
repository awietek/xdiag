#include "tj_symmetric_apply.h"

#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_exchange.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_hopping.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_ising.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetric<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           tJSymmetric<bit_t> const &block_out,
           lila::Vector<coeff_t> &vec_out) {
  using namespace terms::tj_symmetric;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_operator_works_with<coeff_t>(bonds, couplings, block_in.irrep(),
                                            block_out.irrep(),
                                            "tj_symmetric_apply");
  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_hopping_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_ising_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
}

template void Apply<uint16_t, double>(BondList const &, Couplings const &,
                                      tJSymmetric<uint16_t> const &,
                                      lila::Vector<double> const &,
                                      tJSymmetric<uint16_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint32_t, double>(BondList const &, Couplings const &,
                                      tJSymmetric<uint32_t> const &,
                                      lila::Vector<double> const &,
                                      tJSymmetric<uint32_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint64_t, double>(BondList const &, Couplings const &,
                                      tJSymmetric<uint64_t> const &,
                                      lila::Vector<double> const &,
                                      tJSymmetric<uint64_t> const &,
                                      lila::Vector<double> &);

template void Apply<uint16_t, complex>(BondList const &, Couplings const &,
                                       tJSymmetric<uint16_t> const &,
                                       lila::Vector<complex> const &,
                                       tJSymmetric<uint16_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint32_t, complex>(BondList const &, Couplings const &,
                                       tJSymmetric<uint32_t> const &,
                                       lila::Vector<complex> const &,
                                       tJSymmetric<uint32_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint64_t, complex>(BondList const &, Couplings const &,
                                       tJSymmetric<uint64_t> const &,
                                       lila::Vector<complex> const &,
                                       tJSymmetric<uint64_t> const &,
                                       lila::Vector<complex> &);

} // namespace hydra
