#pragma once
#ifdef HYDRA_USE_MPI

#include <mpi.h>
#include <hydra/common.h>

namespace hydra::mpi {

template <class coeff_t> MPI_Datatype datatype();

} // namespace hydra::mpi

#endif
