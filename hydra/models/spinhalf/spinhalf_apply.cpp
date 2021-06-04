#include "spinhalf_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/spinhalf/terms/spinhalf_hopping.h>
#include <hydra/models/spinhalf/terms/spinhalf_u.h>

namespace hydra {

template <class bit_t>
void apply(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           Spinhalf<bit_t> const &block_out, lila::Vector<double> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, double val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };
  
  spinhalfdetail::do_ising(bonds, couplings, block_in, fill);
  spinhalfdetail::do_exchange<bit_t, complex>(bonds, couplings, block_in, fill);
}

template <class bit_t>
void apply(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, lila::Vector<complex> const &vec_in,
           Spinhalf<bit_t> const &block_out, lila::Vector<complex> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  auto fill = [&vec_out, &vec_in](idx_t idx_out, idx_t idx_in, complex val) {
    vec_out(idx_out) += val * vec_in(idx_in);
  };

  spinhalfdetail::do_ising(bonds, couplings, block_in, fill);
  spinhalfdetail::do_exchange<bit_t, complex>(bonds, couplings, block_in, fill);
}

template void apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint16> const &block_in,
                            lila::Vector<double> const &vec_in,
                            Spinhalf<uint16> const &block_out,
                            lila::Vector<double> &vec_out);
template void apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint32> const &block_in,
                            lila::Vector<double> const &vec_in,
                            Spinhalf<uint32> const &block_out,
                            lila::Vector<double> &vec_out);
template void apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint64> const &block_in,
                            lila::Vector<double> const &vec_in,
                            Spinhalf<uint64> const &block_out,
                            lila::Vector<double> &vec_out);

template void apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint16> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            Spinhalf<uint16> const &block_out,
                            lila::Vector<complex> &vec_out);
template void apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint32> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            Spinhalf<uint32> const &block_out,
                            lila::Vector<complex> &vec_out);
template void apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint64> const &block_in,
                            lila::Vector<complex> const &vec_in,
                            Spinhalf<uint64> const &block_out,
                            lila::Vector<complex> &vec_out);

} // namespace hydra
