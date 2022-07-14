#include "print_mpi.h"
#include <hydra/utils/print.h>

namespace hydra::utils {

void PrintPrettyMPI(const char* identifier, Bond const& bond) 
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0) {
    PrintPretty(identifier, bond);
  }
}

void PrintPrettyMPI(const char* identifier, BondList const& bondlist) 
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0) {
    PrintPretty(identifier, bondlist);
  }
}

void PrintPrettyMPI(const char* identifier, Couplings const& couplings) 
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0) {
    PrintPretty(identifier, couplings);
  }
}
  
void PrintPrettyMPI(const char* identifier, Tmatrix const& tmat){
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0) {
    PrintPretty(identifier, tmat);
  }
}

} // namespace hydra::utils
