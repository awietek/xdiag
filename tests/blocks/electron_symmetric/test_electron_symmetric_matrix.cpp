#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.hpp"
#include <xdiag/blocks/electron/electron_matrix.hpp>
#include <xdiag/blocks/electron/electron_apply.hpp>
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/utils/close.hpp>
#include <xdiag/utils/print_macro.hpp>


using namespace xdiag;

void test_electron_symmetric_spectra_no_np(BondList bondlist,
                                           PermutationGroup space_group,
                                           std::vector<Representation> irreps,
                                           std::vector<int64_t> multiplicities) {
  (void) multiplicities;

  int64_t n_sites = space_group.n_sites();

  auto block_total = Electron(n_sites);
  if (block_total.size() < 1000) {

    auto H_total = matrixC(bondlist, block_total, block_total);
    REQUIRE(H_total.is_hermitian(1e-8));

    arma::vec eigs_total;
    arma::eig_sym(eigs_total, H_total);

    std::vector<double> eigs_all;
    for (auto irrep : irreps) {
      auto block_no_np = Electron(n_sites, space_group, irrep);
      auto H_no_np = matrixC(bondlist, block_no_np, block_no_np);

      REQUIRE(H_no_np.is_hermitian(1e-8));

      arma::vec eigs_no_np;
      arma::eig_sym(eigs_no_np, H_no_np);

      std::vector<double> eigs_np_all;
      for (int64_t nup = 0; nup <= n_sites; ++nup) {
        for (int64_t ndn = 0; ndn <= n_sites; ++ndn) {
          auto block_np = Electron(n_sites, nup, ndn, space_group, irrep);
          auto H_np = matrixC(bondlist, block_np, block_np);
          REQUIRE(H_np.is_hermitian(1e-8));

          arma::vec eigs_np;
          arma::eig_sym(eigs_np, H_np);

          for (auto e : eigs_np) {
            eigs_np_all.push_back(e);
            eigs_all.push_back(e);
          }
        }
      }
      std::sort(eigs_np_all.begin(), eigs_np_all.end());
      // XDIAG_PRINT(eigs_no_np);
      // XDIAG_PRINT(arma::vec(eigs_np_all));

      REQUIRE(close(eigs_no_np, arma::vec(eigs_np_all)));
    }
    std::sort(eigs_all.begin(), eigs_all.end());
    REQUIRE(close(eigs_total, arma::vec(eigs_all)));
  }
}

void test_electron_symmetric_spectra(BondList bondlist,
                                     PermutationGroup space_group,
                                     std::vector<Representation> irreps,
                                     std::vector<int64_t> multiplicities) {

  int64_t n_sites = space_group.n_sites();
  assert(irreps.size() == multiplicities.size());

  for (int64_t nup = 0; nup <= n_sites; ++nup) {
    for (int64_t ndn = 0; ndn <= n_sites; ++ndn) {

      // Compute the full spectrum from non-symmetrized block
      auto electron_nosym = Electron(n_sites, nup, ndn);
      if (electron_nosym.size() < 1000) {

        auto H_nosym = matrixC(bondlist, electron_nosym, electron_nosym);
        REQUIRE(H_nosym.is_hermitian(1e-8));
        arma::vec eigs_nosym;
        arma::eig_sym(eigs_nosym, H_nosym);

        std::vector<double> eigs_sym;
        for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int64_t multiplicity = multiplicities[k];

          auto electron = Electron(n_sites, nup, ndn, space_group, irrep);
          // Log.out(
          //     "nup: {}, ndn: {}, k: {}, mult: {}, dim_nosym: {}, dim_sym:"
          //     "{} ",
          //     nup, ndn, k, multiplicity, electron_nosym.size(),
          //     electron.size());

          if (electron.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = matrixC(bondlist, electron, electron);
            REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);

            // REQUIRE(H_sym.is_hermitian(1e-7));
            arma::vec eigs_sym_k;
            arma::eig_sym(eigs_sym_k, H_sym);

            // Check whether results are the same for real blocks
            if (electron.irrep().isreal() && bondlist.isreal()) {
              auto H_sym_real = matrix(bondlist, electron, electron);
              arma::vec eigs_sym_k_real;
              arma::eig_sym(eigs_sym_k_real, H_sym_real);
              REQUIRE(close(eigs_sym_k, eigs_sym_k_real));
            }

            // append all the eigenvalues with multiplicity
            for (auto eig : eigs_sym_k)
              for (int64_t i = 0; i < multiplicity; ++i)
                eigs_sym.push_back(eig);
          }
        }
        std::sort(eigs_sym.begin(), eigs_sym.end());

        // Check if all eigenvalues agree
        // Log.out("{} {} {} {}", nup, ndn, eigs_sym(0), eigs_nosym(0));

	// if (!close(arma::vec(eigs_sym), eigs_nosym)){
	//   XDIAG_PRINT(arma::norm(arma::vec(eigs_sym) - eigs_nosym));
	//   XDIAG_PRINT(eigs_sym);
	//   XDIAG_PRINT(eigs_nosym);
	// }
        REQUIRE(close(arma::vec(eigs_sym), eigs_nosym));
      }
    }
  }
}

void test_hubbard_symmetric_spectrum_chains(int64_t n_sites) {
  using namespace xdiag::testcases::electron;

  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);

  // Without Heisenberg term
  Log.out("electron_symmetric_matrix: Hubbard chain, n_sites: {}", n_sites);
  auto bondlist = get_linear_chain(n_sites, 1.0, 5.0);
  test_electron_symmetric_spectra(bondlist, space_group, irreps,
                                  multiplicities);
  test_electron_symmetric_spectra_no_np(bondlist, space_group, irreps,
                                        multiplicities);

  // With Heisenberg term
  Log("electron_symmetric_matrix: Hubbard chain, n_sites: {} (+ "
      "Heisenberg terms)",
      n_sites);
  auto bondlist_hb = get_linear_chain_hb(n_sites, 0.4);
  test_electron_symmetric_spectra(bondlist_hb, space_group, irreps,
                                  multiplicities);
  test_electron_symmetric_spectra_no_np(bondlist_hb, space_group, irreps,
                                        multiplicities);
}

TEST_CASE("electron_symmetric_matrix", "[electron]") {
  using namespace xdiag::testcases::electron;
  //////////////////////////////////////////////////////////////////////////////////////

  // Check matrices agains Weisse & Fehske
  Log("electron_symmetric_matrix: Weisse & Fehske matrix");
  int64_t n_sites = 4;
  int64_t nup = 3;
  int64_t ndn = 2;
  double t = 1.0;
  double U = 5.0;
  auto bondlist = get_linear_chain(n_sites, t, U);
  bondlist["U"] = U;
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);

  for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
    auto irrep = irreps[k];
    auto electron = Electron(n_sites, nup, ndn, space_group, irrep);
    auto H_sym = matrixC(bondlist, electron, electron);
    complex U2 = 2 * U;
    complex UU = U;
    complex tp = t;
    complex tm = -t;
    complex it = complex(0, t);
    if (k == 0) {
      arma::Mat<complex> H_correct = {
          {U2, tm, tm, tp, tp, 0.}, {tm, U2, tm, tm, 0., tp},
          {tm, tm, U2, 0., tm, tm}, {tp, tm, 0., UU, tm, tp},
          {tp, 0., tm, tm, UU, tm}, {0., tp, tm, tp, tm, UU}};
      REQUIRE(close(H_correct, H_sym));
    }
    if (k == 1) {
      arma::Mat<complex> H_correct = {
          {U2, tm, it, it, tp, 0.},       {tm, U2, tm, tm, 2. * it, tp},
          {-it, tm, U2, 0., tm, it},      {-it, tm, 0., UU, tm, it},
          {tp, -2. * it, tm, tm, UU, tm}, {0., tp, -it, -it, tm, UU}};
      REQUIRE(close(H_correct, H_sym));
    }
    if (k == 2) {
      arma::Mat<complex> H_correct = {
          {U2, tm, tp, tm, tp, 0.}, {tm, U2, tm, tm, 0., tp},
          {tp, tm, U2, 0., tm, tp}, {tm, tm, 0., UU, tm, tm},
          {tp, 0., tm, tm, UU, tm}, {0., tp, tp, tm, tm, UU}};
      REQUIRE(close(H_correct, H_sym));
    }
    if (k == 3) {
      arma::Mat<complex> H_correct = {
          {U2, tm, -it, -it, tp, 0.},    {tm, U2, tm, tm, -2. * it, tp},
          {it, tm, U2, 0., tm, -it},     {it, tm, 0., UU, tm, -it},
          {tp, 2. * it, tm, tm, UU, tm}, {0., tp, it, it, tm, UU}};
      REQUIRE(close(H_correct, H_sym));
    }
  }

  // Test linear chains
  for (int64_t n_sites = 2; n_sites < 7; ++n_sites) {
    test_hubbard_symmetric_spectrum_chains(n_sites);
    test_hubbard_symmetric_spectrum_chains(n_sites);
    test_hubbard_symmetric_spectrum_chains(n_sites);
  }

  // test a 3x3 triangular lattice
  Log("electron_symmetric_matrix: Hubbard 3x3 triangular");
  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.lat";

  bondlist = read_bondlist(lfile);
  bondlist["T"] = 1.0;
  bondlist["U"] = 5.0;
  auto permutations = xdiag::read_permutations(lfile);
  space_group = PermutationGroup(permutations);

  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};
  irreps.clear();
  multiplicities.clear();
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_representation(lfile, name));
    multiplicities.push_back(mult);
  }
  test_electron_symmetric_spectra(bondlist, space_group, irreps,
                                  multiplicities);

  // test a 3x3 triangular lattice with Heisenberg terms
  Log("electron_symmetric_matrix: Hubbard 3x3 triangular(+ Heisenberg terms)");
  auto bondlist_hb = bondlist;
  for (auto bond : bondlist) {
    bondlist_hb << Bond("HB", "J", {bond[0], bond[1]});
  }
  bondlist_hb["J"] = 0.4;
  test_electron_symmetric_spectra(bondlist_hb, space_group, irreps,
                                  multiplicities);

  // test a 3x3 triangular lattice with complex hoppings
  {
    Log.out("electron_symmetric_matrix: Hubbard 3x3 triangular (complex)");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";
    BondList bondlist = read_bondlist(lfile);
    bondlist["TPHI"] = complex(0.5, 0.5);
    bondlist["JPHI"] = 0.;
    bondlist["U"] = 5.0;
    auto permutations = xdiag::read_permutations(lfile);
    space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};
    irreps.clear();
    multiplicities.clear();
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(lfile, name));
      multiplicities.push_back(mult);
    }
    test_electron_symmetric_spectra(bondlist, space_group, irreps,
                                    multiplicities);
  }
}
