#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t>
void apply(BondList const &bonds, Couplings const &couplings,
           SpinhalfMPI<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           SpinhalfMPI<bit_t> const &block_out, lila::Vector<double> &vec_out);

} // namespace hydra
