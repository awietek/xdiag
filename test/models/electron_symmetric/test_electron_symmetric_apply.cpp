#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

template <class bit_t>
void test_symmetric_apply(BondList bondlist, Couplings couplings,
                          SpaceGroup<bit_t> space_group,
                          std::vector<Representation> irreps) {
  int n_sites = space_group.n_sites();

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      for (auto irrep : irreps) {

        // Create block and matrix for comparison
        auto block =
            ElectronSymmetric<bit_t>(n_sites, nup, ndn, space_group, irrep);
        if (block.size() > 0) {
          auto H = MatrixCplx(bondlist, couplings, block, block);
          REQUIRE(lila::close(H, lila::Herm(H)));

          // Check whether apply gives the same as matrix multiplication
          auto v = lila::Random<complex>(block.size());
          auto w1 = lila::Mult(H, v);
          auto w2 = lila::ZerosLike(v);
          Apply(bondlist, couplings, block, v, block, w2);
          REQUIRE(lila::close(w1, w2));

          // Compute eigenvalues and compare
          auto evals_mat = lila::EigenvaluesSym(H);
          double e0_mat = evals_mat(0);
          double e0_app = E0Cplx(bondlist, couplings, block);
          // lila::Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
          REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);
        }
      }
    }
  }
}

template <class bit_t> void test_hubbard_symmetric_apply_chains(int n_sites) {
  using namespace hydra::testcases::electron;
  lila::Log.out("Hubbard chain, symmetric apply test, n_sites: {}", n_sites);
  auto [bondlist, couplings] = get_linear_chain(n_sites, 1.0, 5.0);
  auto [space_group, irreps] = get_cyclic_group_irreps<bit_t>(n_sites);
  test_symmetric_apply(bondlist, couplings, space_group, irreps);
}

TEST_CASE("electron_symmetric_apply", "[models]") {

  // Test linear chains
  for (int n_sites = 2; n_sites < 7; ++n_sites) {
    test_hubbard_symmetric_apply_chains<hydra::uint16>(n_sites);
    test_hubbard_symmetric_apply_chains<hydra::uint32>(n_sites);
    test_hubbard_symmetric_apply_chains<hydra::uint64>(n_sites);
  }

  // test a 3x3 triangular lattice
  lila::Log.out("Hubbard 3x3 triangular, symmetric apply test");
  using bit_t = uint16;
  std::string lfile = "data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";

  auto bondlist = read_bondlist(lfile);
  Couplings couplings;
  couplings["T"] = 1.0;
  couplings["U"] = 5.0;
  auto permutations = read_permutations(lfile);
  auto space_group = SpaceGroup<bit_t>(permutations);

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};
  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_represenation(lfile, name));
  }
  test_symmetric_apply(bondlist, couplings, space_group, irreps);
}
