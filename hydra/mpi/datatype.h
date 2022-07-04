#pragma once

#ifdef HYDRA_ENABLE_MPI

#include <hydra/common.h>
#include <mpi.h>

namespace hydra::mpi {

template <class TCoeffs> MPI_Datatype datatype();

} // namespace hydra::mpi

#endif
