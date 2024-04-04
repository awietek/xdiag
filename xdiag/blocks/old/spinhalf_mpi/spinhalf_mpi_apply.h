#pragma once
#ifdef XDIAG_ENABLE_MPI

#include <mpi.h>
#include <lila/all.h>

#include <xdiag/common.h>
#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/operators/couplings.h>

namespace xdiag {

template <class bit_t, class coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfMPI<bit_t> const &block_in,
           lila::Vector<coeff_t> const &vec_in,
           SpinhalfMPI<bit_t> const &block_out, lila::Vector<coeff_t> &vec_out);

} // namespace xdiag
#endif
