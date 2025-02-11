#include <xdiag/all.hpp>

using namespace xdiag;

int main(int argc, char *argv[]) try {
  assert(argc == 2);
  int64_t nsites = atoi(argv[1]);
  int64_t nup = nsites / 2;

  say_hello();
  set_verbosity(2);

  assert((nsites % 4 == 2) || (nsites % 4 == 0));
  if (nsites % 4 == 0) {
    auto fl =
        FileToml(fmt::format("lattice-files/chain.{}.J1J2.4sl.toml", nsites));
    auto ops = read_opsum(fl, "Interactions");
    auto irrep = read_representation(fl, "Gamma.C2.A");
    ops["J1"] = 1.0;
    ops["J2"] = 0.0;
    tic();
    auto block = Spinhalf(nsites, nup, irrep);
    // auto block = Spinhalf(nsites, nup, irrep, "4sublattice");
    toc("Block creation");

    // tic();
    // double e0 = eigval0(ops, block, 1e-12, 1);
    // toc("MVM");

  } else {
    auto fl =
        FileToml(fmt::format("lattice-files/chain.{}.J1J2.2sl.toml", nsites));
    auto ops = read_opsum(fl, "Interactions");
    auto irrep = read_representation(fl, "X.C2.A");
    ops["J1"] = 1.0;
    ops["J2"] = 0.0;

    tic();
    auto block = Spinhalf(nsites, nup, irrep);
    // auto block = Spinhalf(nsites, nup, irrep, "2sublattice");
    toc("Block creation");

    tic();
    double e0 = eigval0(ops, block, 1e-12, 1);
    toc("MVM");
  }
} catch (Error e) {
  error_trace(e);
}
