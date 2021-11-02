#include "tj_apply.h"

#include <hydra/blocks/tj/terms/tj_exchange.h>
#include <hydra/blocks/tj/terms/tj_hopping.h>
#include <hydra/blocks/tj/terms/tj_ising.h>

namespace hydra {

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJ<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           tJ<bit_t> const &block_out, lila::Vector<double> &vec_out) {
  using namespace terms::tj;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());
  utils::check_operator_real(bonds, couplings, "apply real tJ operator");

  lila::Zeros(vec_out);

  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_hopping<bit_t, double>(bonds, couplings, block_in, fill);
  do_ising<bit_t>(bonds, couplings, block_in, fill);
  do_exchange<bit_t>(bonds, couplings, block_in, fill);
}

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJ<bit_t> const &block_in, lila::Vector<complex> const &vec_in,
           tJ<bit_t> const &block_out, lila::Vector<complex> &vec_out) {
  using namespace terms::tj;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);

  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_hopping<bit_t, complex>(bonds, couplings, block_in, fill);
  do_ising<bit_t>(bonds, couplings, block_in, fill);
  do_exchange<bit_t>(bonds, couplings, block_in, fill);
}

template void Apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint16> const &block_in,
                            lila::Vector<double> const &vec_in,
                            tJ<uint16> const &block_out,
                            lila::Vector<double> &vec_out);
template void Apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint32> const &block_in,
                            lila::Vector<double> const &vec_in,
                            tJ<uint32> const &block_out,
                            lila::Vector<double> &vec_out);
template void Apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint64> const &block_in,
                            lila::Vector<double> const &vec_in,
                            tJ<uint64> const &block_out,
                            lila::Vector<double> &vec_out);

template void Apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint16> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            tJ<uint16> const &block_out,
                            lila::Vector<complex> &vec_out);
template void Apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint32> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            tJ<uint32> const &block_out,
                            lila::Vector<complex> &vec_out);
template void Apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            tJ<uint64> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            tJ<uint64> const &block_out,
                            lila::Vector<complex> &vec_out);

} // namespace hydra
