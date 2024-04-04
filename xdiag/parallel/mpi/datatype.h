#pragma once
#ifdef XDIAG_USE_MPI

#include <mpi.h>
#include <xdiag/common.h>

namespace xdiag::mpi {

template <class coeff_t> MPI_Datatype datatype();

} // namespace xdiag::mpi

#endif
