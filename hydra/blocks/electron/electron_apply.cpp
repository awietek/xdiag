#include "electron_apply.h"

#include <hydra/blocks/electron/terms/electron_exchange.h>
#include <hydra/blocks/electron/terms/electron_hopping.h>
#include <hydra/blocks/electron/terms/electron_ising.h>
#include <hydra/blocks/electron/terms/electron_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           Electron<bit_t> const &block_out, lila::Vector<double> &vec_out) {
  using namespace terms::electron;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  utils::check_operator_real(bonds, couplings, "apply real Electron operator");

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_U(couplings, block_in, fill);
  do_hopping<bit_t, double>(bonds, couplings, block_in, fill);
  do_ising<bit_t>(bonds, couplings, block_in, fill);
  do_exchange<bit_t>(bonds, couplings, block_in, fill);
}

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, lila::Vector<complex> const &vec_in,
           Electron<bit_t> const &block_out, lila::Vector<complex> &vec_out) {
  using namespace terms::electron;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  do_U(couplings, block_in, fill);
  do_hopping<bit_t, complex>(bonds, couplings, block_in, fill);
  do_ising<bit_t>(bonds, couplings, block_in, fill);
  do_exchange<bit_t>(bonds, couplings, block_in, fill);
}

template void Apply<uint16_t>(BondList const &bonds, Couplings const &couplings,
                              Electron<uint16_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              Electron<uint16_t> const &block_out,
                              lila::Vector<double> &vec_out);
template void Apply<uint32_t>(BondList const &bonds, Couplings const &couplings,
                              Electron<uint32_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              Electron<uint32_t> const &block_out,
                              lila::Vector<double> &vec_out);
template void Apply<uint64_t>(BondList const &bonds, Couplings const &couplings,
                              Electron<uint64_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              Electron<uint64_t> const &block_out,
                              lila::Vector<double> &vec_out);

template void Apply<uint16_t>(BondList const &bonds, Couplings const &couplings,
                              Electron<uint16_t> const &block_in,
                              lila::Vector<complex> const &vec_in,
                              Electron<uint16_t> const &block_out,
                              lila::Vector<complex> &vec_out);
template void Apply<uint32_t>(BondList const &bonds, Couplings const &couplings,
                              Electron<uint32_t> const &block_in,
                              lila::Vector<complex> const &vec_in,
                              Electron<uint32_t> const &block_out,
                              lila::Vector<complex> &vec_out);
template void Apply<uint64_t>(BondList const &bonds, Couplings const &couplings,
                              Electron<uint64_t> const &block_in,
                              lila::Vector<complex> const &vec_in,
                              Electron<uint64_t> const &block_out,
                              lila::Vector<complex> &vec_out);

} // namespace hydra
