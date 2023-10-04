#pragma once

#ifdef HYDRA_USE_MPI
#include <mpi.h>
#include <hydra/parallel/mpi/logger_mpi.h>
#else
#include <hydra/utils/logger_serial.h>
#endif

namespace hydra {

#ifdef HYDRA_USE_MPI
inline auto & Log = LogMPI;
#else
inline auto & Log = LogSerial;
#endif

}
