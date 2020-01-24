#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <hydra/utils/combinatorics.h>
#include <lime/all.h>

lila::LoggerMPI lg;

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;
  uint64 downspins = 37;
  uint64 holes = 2;

  std::cout << hydra::combinatorics::down_hole_to_up(downspins, holes) << std::endl;
}
