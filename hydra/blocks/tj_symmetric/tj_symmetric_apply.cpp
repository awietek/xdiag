#include "tj_symmetric_apply.h"

#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_exchange.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_hopping.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_ising.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetric<bit_t> const &block_in,
           lila::Vector<double> const &vec_in,
           tJSymmetric<bit_t> const &block_out, lila::Vector<double> &vec_out) {
  using namespace terms::tj_symmetric;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_symmetric_operator_real(bonds, couplings, block_in.irrep(),
                                       block_out.irrep(),
                                       "apply real tJSymmetric operator");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_hopping_symmetric<bit_t, double>(bonds, couplings, indexing_in, fill);
  do_ising_symmetric<bit_t, double>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, double>(bonds, couplings, indexing_in, fill);
}

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetric<bit_t> const &block_in,
           lila::Vector<complex> const &vec_in,
           tJSymmetric<bit_t> const &block_out,
           lila::Vector<complex> &vec_out) {
  using namespace terms::tj_symmetric;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };
  auto const &indexing_in = block_in.indexing();

  do_hopping_symmetric<bit_t, complex>(bonds, couplings, indexing_in, fill);
  do_ising_symmetric<bit_t, complex>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, complex>(bonds, couplings, indexing_in, fill);
}

template void Apply<uint16_t>(BondList const &, Couplings const &,
                              tJSymmetric<uint16_t> const &,
                              lila::Vector<double> const &,
                              tJSymmetric<uint16_t> const &,
                              lila::Vector<double> &);
template void Apply<uint32_t>(BondList const &, Couplings const &,
                              tJSymmetric<uint32> const &,
                              lila::Vector<double> const &,
                              tJSymmetric<uint32> const &,
                              lila::Vector<double> &);
template void Apply<uint64_t>(BondList const &, Couplings const &,
                              tJSymmetric<uint64> const &,
                              lila::Vector<double> const &,
                              tJSymmetric<uint64> const &,
                              lila::Vector<double> &);

template void Apply<uint16_t>(BondList const &, Couplings const &,
                              tJSymmetric<uint16_t> const &,
                              lila::Vector<complex> const &,
                              tJSymmetric<uint16_t> const &,
                              lila::Vector<complex> &);
template void Apply<uint32_t>(BondList const &, Couplings const &,
                              tJSymmetric<uint32> const &,
                              lila::Vector<complex> const &,
                              tJSymmetric<uint32> const &,
                              lila::Vector<complex> &);
template void Apply<uint64_t>(BondList const &, Couplings const &,
                              tJSymmetric<uint64> const &,
                              lila::Vector<complex> const &,
                              tJSymmetric<uint64> const &,
                              lila::Vector<complex> &);

} // namespace hydra
