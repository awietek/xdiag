#include "spinhalf_mpi_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/spinhalf_mpi/terms/spinhalf_mpi_ising.h>
// #include <hydra/models/spinhalf_mpi/terms/spinhalf_mpi_exchange.h>

namespace hydra {

template <class bit_t>
void apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfMPI<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           SpinhalfMPI<bit_t> const &block_out, lila::Vector<double> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);
  
  spinhalfterms::do_ising_mpi(bonds, couplings, block_in, vec_in, vec_out);
  // spinhalfterms::do_exchange_mpi(bonds, couplings, block_in);
}

template void apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            SpinhalfMPI<uint16> const &block_in,
                            lila::Vector<double> const &vec_in,
                            SpinhalfMPI<uint16> const &block_out,
                            lila::Vector<double> &vec_out);
template void apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            SpinhalfMPI<uint32> const &block_in,
                            lila::Vector<double> const &vec_in,
                            SpinhalfMPI<uint32> const &block_out,
                            lila::Vector<double> &vec_out);
template void apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            SpinhalfMPI<uint64> const &block_in,
                            lila::Vector<double> const &vec_in,
                            SpinhalfMPI<uint64> const &block_out,
                            lila::Vector<double> &vec_out);


} // namespace hydra
