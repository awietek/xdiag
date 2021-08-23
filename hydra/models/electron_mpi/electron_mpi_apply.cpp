#include "electron_mpi_apply.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/utils/bitops.h>

#include <hydra/models/electron_mpi/terms/electron_mpi_hopping.h>
#include <hydra/models/electron_mpi/terms/electron_mpi_u.h>

namespace hydra {

template <class bit_t, class coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronMPI<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           ElectronMPI<bit_t> const &block_out,
           lila::Vector<coeff_t> &vec_out) {

  assert(block_in == block_out); // only temporary
  assert(block_in.size() == vec_in.size());
  assert(block_out.size() == vec_out.size());

  lila::Zeros(vec_out);

  auto tis = rightnow_mpi();
  electron::do_U_mpi(bonds, couplings, block_in, vec_in, vec_out);
  timing_mpi(tis, rightnow_mpi(), " (U)", 2);

  // auto tex = rightnow_mpi();
  // electron::do_hopping_mpi(bonds, couplings, block_in, vec_in, vec_out);
  // timing_mpi(tex, rightnow_mpi(), " (exchange)", 2);
}

template void Apply<uint16>(BondList const &bonds, Couplings const &couplings,
                            ElectronMPI<uint16> const &block_in,
                            lila::Vector<double> const &vec_in,
                            ElectronMPI<uint16> const &block_out,
                            lila::Vector<double> &vec_out);
template void Apply<uint32>(BondList const &bonds, Couplings const &couplings,
                            ElectronMPI<uint32> const &block_in,
                            lila::Vector<double> const &vec_in,
                            ElectronMPI<uint32> const &block_out,
                            lila::Vector<double> &vec_out);
template void Apply<uint64>(BondList const &bonds, Couplings const &couplings,
                            ElectronMPI<uint64> const &block_in,
                            lila::Vector<double> const &vec_in,
                            ElectronMPI<uint64> const &block_out,
                            lila::Vector<double> &vec_out);

template void Apply<uint16, complex>(BondList const &bonds,
                                     Couplings const &couplings,
                                     ElectronMPI<uint16> const &block_in,
                                     lila::Vector<complex> const &vec_in,
                                     ElectronMPI<uint16> const &block_out,
                                     lila::Vector<complex> &vec_out);
template void Apply<uint32, complex>(BondList const &bonds,
                                     Couplings const &couplings,
                                     ElectronMPI<uint32> const &block_in,
                                     lila::Vector<complex> const &vec_in,
                                     ElectronMPI<uint32> const &block_out,
                                     lila::Vector<complex> &vec_out);
template void Apply<uint64, complex>(BondList const &bonds,
                                     Couplings const &couplings,
                                     ElectronMPI<uint64> const &block_in,
                                     lila::Vector<complex> const &vec_in,
                                     ElectronMPI<uint64> const &block_out,
                                     lila::Vector<complex> &vec_out);

} // namespace hydra
