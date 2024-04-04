#pragma once
#ifdef XDIAG_ENABLE_MPI

#include <mpi.h>
#include <lila/all.hpp>

#include <xdiag/common.hpp>
#include <xdiag/blocks/electron_mpi/electron_mpi.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/operators/couplings.hpp>

namespace xdiag {

template <class bit_t, class coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronMPI<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           ElectronMPI<bit_t> const &block_out, lila::Vector<coeff_t> &vec_out);

} // namespace xdiag
#endif
