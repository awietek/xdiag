#include "spinhalf_symmetric_apply.h"

#include <hydra/blocks/spinhalf_symmetric/spinhalf_symmetric.h>
#include <hydra/blocks/spinhalf_symmetric/terms/spinhalf_symmetric_exchange.h>
#include <hydra/blocks/spinhalf_symmetric/terms/spinhalf_symmetric_ising.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfSymmetric<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           SpinhalfSymmetric<bit_t> const &block_out,
           lila::Vector<coeff_t> &vec_out) {
  using namespace terms::spinhalf_symmetric;
  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_operator_works_with<coeff_t>(bonds, couplings, block_in.irrep(),
                                            block_out.irrep(),
                                            "spinhalf_symmetric_apply");
  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, coeff_t val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_ising_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
}

template void Apply<uint16_t, double>(BondList const &, Couplings const &,
                                      SpinhalfSymmetric<uint16_t> const &,
                                      lila::Vector<double> const &,
                                      SpinhalfSymmetric<uint16_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint32_t, double>(BondList const &, Couplings const &,
                                      SpinhalfSymmetric<uint32_t> const &,
                                      lila::Vector<double> const &,
                                      SpinhalfSymmetric<uint32_t> const &,
                                      lila::Vector<double> &);
template void Apply<uint64_t, double>(BondList const &, Couplings const &,
                                      SpinhalfSymmetric<uint64_t> const &,
                                      lila::Vector<double> const &,
                                      SpinhalfSymmetric<uint64_t> const &,
                                      lila::Vector<double> &);

template void Apply<uint16_t, complex>(BondList const &, Couplings const &,
                                       SpinhalfSymmetric<uint16_t> const &,
                                       lila::Vector<complex> const &,
                                       SpinhalfSymmetric<uint16_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint32_t, complex>(BondList const &, Couplings const &,
                                       SpinhalfSymmetric<uint32_t> const &,
                                       lila::Vector<complex> const &,
                                       SpinhalfSymmetric<uint32_t> const &,
                                       lila::Vector<complex> &);
template void Apply<uint64_t, complex>(BondList const &, Couplings const &,
                                       SpinhalfSymmetric<uint64_t> const &,
                                       lila::Vector<complex> const &,
                                       SpinhalfSymmetric<uint64_t> const &,
                                       lila::Vector<complex> &);

} // namespace hydra
