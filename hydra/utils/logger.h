#pragma once

#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#include <hydra/utils/logger_mpi.h>
#else
#include <hydra/utils/logger_serial.h>
#endif

namespace hydra {

#ifdef HYDRA_ENABLE_MPI
inline auto & Log = LogMPI;
#else
inline auto & Log = LogSerial;
#endif

}
