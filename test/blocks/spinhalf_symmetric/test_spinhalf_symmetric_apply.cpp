#include "../../catch.hpp"

#include "../electron/testcases_electron.h"
#include "../spinhalf/testcases_spinhalf.h"

#include <iostream>

#include <hydra/all.h>
using namespace hydra;

template <class bit_t>
void test_spinhalf_symmetric_apply(BondList bondlist, Couplings couplings,
                                   PermutationGroup space_group,
                                   std::vector<Representation> irreps) {
  int n_sites = space_group.n_sites();

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (auto irrep : irreps) {
      auto block = Spinhalf<bit_t>(n_sites, nup, space_group, irrep);

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

        // Compute eigenvalues with real arithmitic
        if (is_real(block.irrep()) && is_real(couplings)) {
          auto H_real = MatrixReal(bondlist, couplings, block, block);
          auto evals_mat_real = lila::EigenvaluesSym(H_real);
          REQUIRE(lila::close(evals_mat_real, evals_mat));

          double e0_mat_real = evals_mat_real(0);
          double e0_app_real = E0Real(bondlist, couplings, block);
          REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
        }
      }
    }
  }
}

template <class bit_t>
void test_spinhalf_symmetric_apply_no_sz(BondList bondlist, Couplings couplings,
                                         PermutationGroup space_group,
                                         std::vector<Representation> irreps) {
  int n_sites = space_group.n_sites();

  for (auto irrep : irreps) {
    auto block = Spinhalf<bit_t>(n_sites, space_group, irrep);

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

      // Compute eigenvalues with real arithmitic
      if (is_real(block.irrep()) && is_real(couplings)) {
        auto H_real = MatrixReal(bondlist, couplings, block, block);
        auto evals_mat_real = lila::EigenvaluesSym(H_real);
        REQUIRE(lila::close(evals_mat_real, evals_mat));

        double e0_mat_real = evals_mat_real(0);
        double e0_app_real = E0Real(bondlist, couplings, block);
        REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
      }
    }
  }
}

template <class bit_t> void test_spinhalf_symmetric_apply_chains(int n_sites) {
  using namespace hydra::testcases::spinhalf;
  using hydra::testcases::electron::get_cyclic_group_irreps;
  lila::Log.out("spinhalf_symmetric_apply: HB chain, N: {}", n_sites);
  auto [space_group, irreps] = get_cyclic_group_irreps(n_sites);
  auto [bondlist, couplings] = HBchain(n_sites, 1.0, 1.0);
  test_spinhalf_symmetric_apply<bit_t>(bondlist, couplings, space_group,
                                       irreps);
  test_spinhalf_symmetric_apply_no_sz<bit_t>(bondlist, couplings, space_group,
                                             irreps);
}

TEST_CASE("spinhalf_symmetric_apply", "[blocks][spinhalf_symmetric]") {

  // Test linear Heisenberg chains
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    test_spinhalf_symmetric_apply_chains<uint16_t>(n_sites);
    test_spinhalf_symmetric_apply_chains<uint32_t>(n_sites);
    test_spinhalf_symmetric_apply_chains<uint64_t>(n_sites);
  }

  // test a 3x3 triangular lattice
  lila::Log.out("spinhalf_symmetric_apply: Triangular 3x3");
  std::string lfile = "data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";

  auto bondlist = read_bondlist(lfile);
  Couplings couplings;
  couplings["Jz1"] = 1.00;
  couplings["Jz2"] = 0.23;
  couplings["Jx1"] = 0.76;
  couplings["Jx2"] = 0.46;

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
      {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
      {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
      {"Y.D1.A", 6},      {"Y.D1.B", 6}};

  auto permutations = hydra::read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);

  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_represenation(lfile, name));
  }
  test_spinhalf_symmetric_apply<uint16_t>(bondlist, couplings, space_group,
                                          irreps);
  test_spinhalf_symmetric_apply_no_sz<uint64_t>(bondlist, couplings,
                                                space_group, irreps);
}
