#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

template <class bit_t>
void test_electron_symmetric_apply(BondList bondlist,
                                   PermutationGroup space_group,
                                   std::vector<Representation> irreps) {
  int n_sites = space_group.n_sites();

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      for (auto irrep : irreps) {

        // Create block and matrix for comparison
        auto block = Electron<bit_t>(n_sites, nup, ndn, space_group, irrep);
        if (block.size() > 0) {
          auto H = matrix_cplx(bondlist, block, block);
          REQUIRE(arma::norm(H - H.t()) < 1e-12);

          // REQUIRE(H.is_hermitian(1e-8));
          // Check whether apply gives the same as matrix multiplication
          arma::cx_vec v(block.size(), arma::fill::randn);
          arma::cx_vec w1 = H * v;
          arma::cx_vec w2(block.size(), arma::fill::randn);
          apply(bondlist, block, v, block, w2);
          REQUIRE(close(w1, w2));
          // Compute eigenvalues and compare
          arma::vec evals_mat;
          arma::eig_sym(evals_mat, H);
          double e0_mat = evals_mat(0);
          double e0_app = e0_cplx(bondlist, block);
          // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
          REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);

          // Compute eigenvalues with real arithmitic
          if (is_real(block.irrep()) && bondlist.is_real()) {
            auto H_real = matrix_real(bondlist, block, block);
            arma::vec evals_mat_real;
            arma::eig_sym(evals_mat_real, H_real);

            REQUIRE(close(evals_mat_real, evals_mat));

            double e0_mat_real = evals_mat_real(0);
            double e0_app_real = e0_real(bondlist, block);
            REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
          }
        }
      }
    }
  }
}

template <class bit_t> void test_hubbard_symmetric_apply_chains(int n_sites) {
  using namespace hydra::testcases::electron;

  // Without Heisenberg term
  Log.out("electron_symmetric_apply: Hubbard chain, n_sites: {}", n_sites);
  auto bondlist = get_linear_chain(n_sites, 1.0, 5.0);
  auto [space_group, irreps] = get_cyclic_group_irreps(n_sites);
  test_electron_symmetric_apply<uint16_t>(bondlist, space_group, irreps);

  // With Heisenberg term
  Log.out("electron_symmetric_apply: Hubbard chain, n_sites: {} (+ "
          "Heisenberg terms)",
          n_sites);
  auto bondlist_hb = get_linear_chain_hb(n_sites, 1.0, 5.0, 0.4);
  test_electron_symmetric_apply<bit_t>(bondlist_hb, space_group, irreps);
}

TEST_CASE("electron_symmetric_apply", "[blocks][electron_symmetric]") {

  // Test linear chains
  for (int n_sites = 2; n_sites < 7; ++n_sites) {
    test_hubbard_symmetric_apply_chains<uint16_t>(n_sites);
    test_hubbard_symmetric_apply_chains<uint32_t>(n_sites);
    test_hubbard_symmetric_apply_chains<uint64_t>(n_sites);
  }

  // test a 3x3 triangular lattice
  Log.out("electron_symmetric_apply: Hubbard 3x3 triangular");
  using bit_t = uint16_t;
  std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";

  auto bondlist = read_bondlist(lfile);
  bondlist["T"] = 1.0;
  bondlist["U"] = 5.0;
  auto permutations = hydra::read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};
  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_represenation(lfile, name));
    (void)mult;
  }
  test_electron_symmetric_apply<bit_t>(bondlist, space_group, irreps);

  // test a 3x3 triangular lattice with Heisenberg terms
  Log.out(
      "electron_symmetric_apply: Hubbard 3x3 triangular (+ Heisenberg terms)");
  auto bondlist_hb = bondlist;
  for (auto bond : bondlist) {
    bondlist_hb << Bond("HB", "J", {bond[0], bond[1]});
  }
  bondlist_hb["J"] = 0.4;
  test_electron_symmetric_apply<bit_t>(bondlist_hb, space_group, irreps);

  // test a 3x3 triangular lattice with complex hoppings
  {
    Log.out("electron_symmetric_apply: Hubbard 3x3 triangular (complex)");
    using bit_t = uint16_t;
    std::string lfile =
        "data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";
    BondList bondlist = read_bondlist(lfile);
    bondlist["TPHI"] = complex(0.5, 0.5);
    bondlist["JPHI"] = 0.;
    bondlist["U"] = 5.0;
    auto permutations = hydra::read_permutations(lfile);
    space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};
    irreps.clear();
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_represenation(lfile, name));
      (void)mult;
    }
    test_electron_symmetric_apply<bit_t>(bondlist, space_group, irreps);
  }
}
