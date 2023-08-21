#pragma once
#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#include <hydra/mpi/logger_mpi.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/linalg/lanczos/tmatrix.h>

namespace hydra::utils {

  void print_pretty_mpi(const char* identifier, Bond const& bond);
  void print_pretty_mpi(const char* identifier, BondList const& bondlist);
  void print_pretty_mpi(const char* identifier, Couplings const& couplings);

  template <typename T>
  void print_pretty_mpi(const char* identifier, lila::Vector<T> const& vector){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0) {
      PrintPretty(identifier, vector);
    }
  }
  template <typename T>
  void print_pretty_mpi(const char* identifier, lila::Matrix<T> const& matrix){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0) {
      PrintPretty(identifier, matrix);
    }
  }

  void print_pretty_mpi(const char* identifier, Tmatrix const& tmat);
  
} // namespace hydra::utils
#endif
