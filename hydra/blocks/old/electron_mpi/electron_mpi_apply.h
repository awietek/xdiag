#pragma once
#ifdef HYDRA_ENABLE_MPI

#include <mpi.h>
#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/blocks/electron_mpi/electron_mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronMPI<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           ElectronMPI<bit_t> const &block_out, lila::Vector<coeff_t> &vec_out);

} // namespace hydra
#endif
