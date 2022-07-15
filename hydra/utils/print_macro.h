#pragma once

#ifdef HYDRA_ENABLE_MPI

#include <mpi.h>
#include <hydra/utils/print_mpi.h>
#define HydraPrint(X) hydra::utils::PrintPrettyMPI(#X,X)

#else

#include <hydra/utils/print.h>
#define HydraPrint(X) hydra::utils::PrintPretty(#X,X)

#endif
