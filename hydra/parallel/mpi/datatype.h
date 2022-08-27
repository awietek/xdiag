#pragma once
#ifdef HYDRA_ENABLE_MPI

#include <mpi.h>
#include <hydra/common.h>

namespace hydra::mpi {

template <class TCoeffs> MPI_Datatype datatype();

} // namespace hydra::mpi

#endif
