#include "electron_matrix.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron/terms/electron_hopping.h>
#include <hydra/models/electron/terms/electron_u.h>

namespace hydra {

template <class bit_t>
void apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           Electron<bit_t> const &block_out, lila::Vector<double> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);

  // Hubbard U
  electron::do_U(couplings, block_in,
                 [&vec_out, &vec_in](idx_t idx, double val) {
                   vec_out(idx) += val * vec_in(idx);
                 });

  // electron hopping
  electron::do_hopping<bit_t, double>(
      bonds, couplings, block_in,
      [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
        vec_out(idx_out) += val * vec_in(idx_in);
	assert(idx_out < vec_out.size());
      });
}

template <class bit_t>
void apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, lila::Vector<complex> const &vec_in,
           Electron<bit_t> const &block_out, lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);

  // Hubbard U
  electron::do_U(couplings, block_in,
                 [&vec_out, &vec_in](idx_t idx, double val) {
                   vec_out(idx) += val * vec_in(idx);
                 });

  // electron hopping
  electron::do_hopping<bit_t, complex>(
      bonds, couplings, block_in,
      [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
        vec_out(idx_out) += val * vec_in(idx_in);
      });
}

template void apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint16> const &block_in,
                            lila::Vector<double> const &vec_in,
                            Electron<uint16> const &block_out,
                            lila::Vector<double> &vec_out);
template void apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint32> const &block_in,
                            lila::Vector<double> const &vec_in,
                            Electron<uint32> const &block_out,
                            lila::Vector<double> &vec_out);
template void apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint64> const &block_in,
                            lila::Vector<double> const &vec_in,
                            Electron<uint64> const &block_out,
                            lila::Vector<double> &vec_out);

template void apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint16> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            Electron<uint16> const &block_out,
                            lila::Vector<complex> &vec_out);
template void apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint32> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            Electron<uint32> const &block_out,
                            lila::Vector<complex> &vec_out);
template void apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint64> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            Electron<uint64> const &block_out,
                            lila::Vector<complex> &vec_out);

} // namespace hydra
