#include "spinhalf_mpi_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/bitops/bitops.h>

#include <hydra/blocks/spinhalf_mpi/terms/spinhalf_mpi_exchange.h>
#include <hydra/blocks/spinhalf_mpi/terms/spinhalf_mpi_ising.h>

namespace hydra {

template <class bit_t, class coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfMPI<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           SpinhalfMPI<bit_t> const &block_out,
           lila::Vector<coeff_t> &vec_out) {
  using namespace terms::spinhalf_mpi;

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);

  auto tis = rightnow_mpi();
  do_ising_mpi(bonds, couplings, block_in, vec_in, vec_out);
  timing_mpi(tis, rightnow_mpi(), " (ising)", 2);

  auto tex = rightnow_mpi();
  do_exchange_mpi(bonds, couplings, block_in, vec_in, vec_out);
  timing_mpi(tex, rightnow_mpi(), " (exchange)", 2);
}

template void Apply<uint16_t>(BondList const &bonds, Couplings const &couplings,
                              SpinhalfMPI<uint16_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              SpinhalfMPI<uint16_t> const &block_out,
                              lila::Vector<double> &vec_out);
template void Apply<uint32_t>(BondList const &bonds, Couplings const &couplings,
                              SpinhalfMPI<uint32_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              SpinhalfMPI<uint32_t> const &block_out,
                              lila::Vector<double> &vec_out);
template void Apply<uint64_t>(BondList const &bonds, Couplings const &couplings,
                              SpinhalfMPI<uint64_t> const &block_in,
                              lila::Vector<double> const &vec_in,
                              SpinhalfMPI<uint64_t> const &block_out,
                              lila::Vector<double> &vec_out);

template void Apply<uint16_t, complex>(BondList const &bonds,
                                       Couplings const &couplings,
                                       SpinhalfMPI<uint16_t> const &block_in,
                                       lila::Vector<complex> const &vec_in,
                                       SpinhalfMPI<uint16_t> const &block_out,
                                       lila::Vector<complex> &vec_out);
template void Apply<uint32_t, complex>(BondList const &bonds,
                                       Couplings const &couplings,
                                       SpinhalfMPI<uint32_t> const &block_in,
                                       lila::Vector<complex> const &vec_in,
                                       SpinhalfMPI<uint32_t> const &block_out,
                                       lila::Vector<complex> &vec_out);
template void Apply<uint64_t, complex>(BondList const &bonds,
                                       Couplings const &couplings,
                                       SpinhalfMPI<uint64_t> const &block_in,
                                       lila::Vector<complex> const &vec_in,
                                       SpinhalfMPI<uint64_t> const &block_out,
                                       lila::Vector<complex> &vec_out);

} // namespace hydra