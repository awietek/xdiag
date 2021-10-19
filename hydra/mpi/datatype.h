#pragma once

#include <hydra/common.h>
#include <mpi.h>

namespace hydra::mpi {

template <class TCoeffs> MPI_Datatype datatype();

} // namespace hydra::mpi
