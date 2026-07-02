// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"
#include <tests/is_approx_hermitian.hpp>

#include <tests/blocks/random_opsum_matrix.hpp>

#include "test_spinhalf_strategies.hpp"
#include "testcases_spinhalf.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/binomial.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/kernels/sparse/csc_matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;
using namespace xdiag::testcases::spinhalf;

// -----------------------------------------------------------------------
// TEST_CASE: basic Spinhalf matrix and sparse matrix construction
// -----------------------------------------------------------------------

TEST_CASE("spinhalf", "[spinhalf]") try {
  int64_t nsites = 6;
  auto block = Spinhalf(nsites, nsites / 2);
  auto ops = HBchain(nsites, 1.0);

  auto H = matrix(ops, block);
  arma::vec eigs;
  arma::eig_sym(eigs, H);
  double e0 = eigval0(ops, block);
  REQUIRE(isapprox(e0, eigs[0]));

  auto H_coo = coo_matrix(ops, block);
  auto H_coo_dense = to_dense(H_coo);
  arma::vec eigs_coo;
  arma::eig_sym(eigs_coo, H_coo_dense);
  REQUIRE(isapprox(H_coo_dense, H));
  REQUIRE(isapprox(e0, eigs_coo[0]));

  auto H_csr = csr_matrix(ops, block);
  auto H_csr_dense = to_dense(H_csr);
  arma::vec eigs_csr;
  arma::eig_sym(eigs_csr, H_csr_dense);
  REQUIRE(isapprox(H_csr_dense, H));
  REQUIRE(isapprox(e0, eigs_csr[0]));

  auto H_csc = csc_matrix(ops, block);
  auto H_csc_dense = to_dense(H_csc);
  arma::vec eigs_csc;
  arma::eig_sym(eigs_csc, H_csc_dense);
  REQUIRE(isapprox(H_csc_dense, H));
  REQUIRE(isapprox(e0, eigs_csc[0]));

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: matrix-vector and matrix-matrix apply
// -----------------------------------------------------------------------

TEST_CASE("spinhalf_apply", "[spinhalf]") try {

  Log("spinhalf_apply: Heisenberg chain apply test, N=2..6");
  for (int64_t N = 2; N <= 6; ++N) {
    auto ops = HBchain(N, 1.0);
    test_apply(N, ops);
    test_apply_mat(N, ops);
  }

  Log("spinhalf_apply: Heisenberg all-to-all apply test, N=2..6");
  for (int64_t N = 2; N <= 6; ++N) {
    auto ops = HB_alltoall(N);
    test_apply(N, ops);
    test_apply_mat(N, ops);
  }

  Log("spinhalf_apply: Heisenberg all-to-all Sz vs. no-Sz ground state");
  for (int64_t nsites = 2; nsites <= 6; ++nsites) {
    auto ops = HB_alltoall(nsites);
    auto block_no_sz = Spinhalf(nsites);
    double e0_no_sz = eigval0(ops, block_no_sz);

    std::vector<double> e0s_sz;
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      auto block_sz = Spinhalf(nsites, nup);
      e0s_sz.push_back(eigval0(ops, block_sz));
    }
    double e0_sz = *std::min_element(e0s_sz.begin(), e0s_sz.end());
    REQUIRE(isapprox(e0_sz, e0_no_sz));
  }

  {
    Log("spinhalf_apply: Triangular J1J2Jchi N=12");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = 0.09;
    int64_t nsites = 12;
    int64_t nup = 6;
    auto spinhalf = Spinhalf(nsites, nup);
    double e0 = eigval0(ops, spinhalf);
    double energy = -6.9456000700824329641;
    REQUIRE(isapprox(e0, energy, 1e-12, 1e-8));
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: dense matrix construction and spectrum verification
// -----------------------------------------------------------------------

TEST_CASE("spinhalf_matrix", "[spinhalf]") try {

  {
    Log("spinhalf_matrix: Heisenberg chain, J=1, N=2..6 exact spectra");
    for (int64_t nsites = 2; nsites <= 6; ++nsites) {
      for (int64_t nup = 0; nup <= nsites; ++nup) {
        auto [ops, exact_eigs] = HBchain_fullspectrum_nup(nsites, nup);
        auto block = Spinhalf(nsites, nup);
        auto H = matrix(ops, block, block);
        REQUIRE(testcases::is_approx_hermitian(H, 1e-7));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        REQUIRE(isapprox(eigs, exact_eigs));
      }
    }
  }

  {
    Log("spinhalf_matrix: Heisenberg all-to-all Sz vs. no-Sz spectrum");
    for (int64_t nsites = 2; nsites <= 6; ++nsites) {
      auto ops = HB_alltoall(nsites);
      auto block_no_sz = Spinhalf(nsites);
      auto H_no_sz = matrix(ops, block_no_sz, block_no_sz);
      REQUIRE(testcases::is_approx_hermitian(H_no_sz, 1e-8));
      arma::vec eigs_no_sz;
      arma::eig_sym(eigs_no_sz, H_no_sz);

      std::vector<double> eigs_sz_all;
      for (int64_t nup = 0; nup <= nsites; ++nup) {
        auto block_sz = Spinhalf(nsites, nup);
        auto H_sz = matrix(ops, block_sz, block_sz);
        REQUIRE(testcases::is_approx_hermitian(H_sz, 1e-7));
        arma::vec eigs_sz;
        arma::eig_sym(eigs_sz, H_sz);
        for (auto eig : eigs_sz) {
          eigs_sz_all.push_back(eig);
        }
      }
      std::sort(eigs_sz_all.begin(), eigs_sz_all.end());
      REQUIRE(isapprox(eigs_no_sz, arma::vec(eigs_sz_all)));
    }
  }

  {
    Log("spinhalf_matrix: triangular N=12 complex exchange");
    int64_t nsites = 12;
    int64_t nup = 6;
    std::vector<double> etas = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05};
    for (auto eta : etas) {
      auto [ops, e0] = triangular_12_complex(nup, eta);
      auto block = Spinhalf(nsites, nup);
      auto H = matrixC(ops, block, block);
      REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
      arma::vec eigs;
      arma::eig_sym(eigs, H);
      REQUIRE(isapprox(eigs(0), e0, 1e-12, 1e-8));
    }
  }

  {
    Log("spinhalf_matrix: Triangular J1J2Jchi N=12");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = 0.09;
    int64_t nsites = 12;
    int64_t nup = 6;
    auto block = Spinhalf(nsites, nup);
    auto H = matrixC(ops, block, block);
    REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
    arma::vec eigs;
    arma::eig_sym(eigs, H);
    double energy = -6.9456000700824329641;
    REQUIRE(isapprox(eigs(0), energy, 1e-12, 1e-8));
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: spin operator algebra and commutation relations
// -----------------------------------------------------------------------

TEST_CASE("spinhalf_commutators", "[spinhalf]") try {

  Log("spinhalf_commutators: [S+_i, S-_j] = 2*Sz_i * delta_ij");
  for (int64_t nsites = 2; nsites < 5; ++nsites) {
    auto block_raw = Spinhalf(nsites);
    for (int64_t nup = 1; nup < nsites; ++nup) {
      auto block = Spinhalf(nsites, nup);
      auto blockp = Spinhalf(nsites, nup + 1);
      auto blockm = Spinhalf(nsites, nup - 1);

      for (int64_t i = 0; i < nsites; ++i) {
        for (int64_t j = 0; j < nsites; ++j) {
          OpSum sp_i_m;
          sp_i_m += Op("S+", i);
          auto sp_i_m_mat = matrix(sp_i_m, blockm, block);
          auto sp_i_mat = matrix(sp_i_m, block_raw, block_raw);

          OpSum sm_j_m;
          sm_j_m += Op("S-", j);
          auto sm_j_m_mat = matrix(sm_j_m, block, blockm);
          auto sm_j_mat = matrix(sm_j_m, block_raw, block_raw);

          OpSum sp_i_p;
          sp_i_p += Op("S+", i);
          auto sp_i_p_mat = matrix(sp_i_p, block, blockp);

          OpSum sm_j_p;
          sm_j_p += Op("S-", j);
          auto sm_j_p_mat = matrix(sm_j_p, blockp, block);

          arma::mat comm = arma::mat(sp_i_m_mat * sm_j_m_mat) -
                           arma::mat(sm_j_p_mat * sp_i_p_mat);
          arma::mat commr =
              arma::mat(sp_i_mat * sm_j_mat) - arma::mat(sm_j_mat * sp_i_mat);

          if (i == j) {
            OpSum sz;
            sz += Op("Sz", i);
            auto sz_mat = matrix(sz, block, block);
            auto sz_matr = matrix(sz, block_raw, block_raw);
            REQUIRE(isapprox(comm, arma::mat(2.0 * sz_mat)));
            REQUIRE(isapprox(commr, arma::mat(2.0 * sz_matr)));
          } else {
            REQUIRE(isapprox(
                comm, arma::mat(comm.n_rows, comm.n_cols, arma::fill::zeros)));
            REQUIRE(isapprox(commr, arma::mat(commr.n_rows, commr.n_cols,
                                              arma::fill::zeros)));
          }
        }
      }
    }
  }

  Log("spinhalf_commutators: onsite Sz*Sz == SzSz");
  test_onsite("Sz", "SzSz");

  {
    Log("spinhalf_commutators: onsite Exchange and SdotS identities");
    for (int64_t nsites = 2; nsites < 5; ++nsites) {
      auto b = Spinhalf(nsites);
      for (int64_t s = 0; s < nsites; ++s) {
        arma::mat sp = matrix(Op("S+", s), b);
        arma::mat sm = matrix(Op("S-", s), b);
        arma::mat sz = matrix(Op("Sz", s), b);

        // Exchange_{s,s} = (1/2)(S+_s S-_s + S-_s S+_s)
        arma::mat m1 = arma::mat(0.5 * (sp * sm + sm * sp));
        arma::mat m2 = matrix(Op("Exchange", {s, s}), b);
        REQUIRE(isapprox(m1, m2));

        // SdotS_{s,s} = Exchange_{s,s} + Sz_s * Sz_s
        m1 = arma::mat(0.5 * (sp * sm + sm * sp) + sz * sz);
        m2 = matrix(Op("SdotS", {s, s}), b);
        REQUIRE(isapprox(m1, m2));
      }
    }
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: symmetry-resolved block dimensions
// -----------------------------------------------------------------------

TEST_CASE("spinhalf_symmetric", "[spinhalf]") try {

  // Cyclic chain: sum of all irrep dimensions equals C(nsites, nup)
  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    Log("spinhalf_symmetric: chain block dimensions, N={}", nsites);
    int64_t sum_of_dims = 0;
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      int64_t sum_of_dims_up = 0;
      for (int64_t k = 0; k < nsites; ++k) {
        auto irrep = cyclic_group_irrep(nsites, k);
        auto block = Spinhalf(nsites, nup, irrep);
        sum_of_dims += block.size();
        sum_of_dims_up += block.size();
      }
      REQUIRE(sum_of_dims_up == math::binomial(nsites, nup));
    }
    REQUIRE(sum_of_dims == ((int64_t)1 << nsites));
  }

  // Triangular 3x3: verify sum of multiplicity-weighted dimensions
  {
    Log("spinhalf_symmetric: triangular 3x3 block dimensions");
    int64_t nsites = 9;
    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
        {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
        {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
        {"Y.D1.A", 6},      {"Y.D1.B", 6}};

    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);

    int64_t sum_dim = 0;
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      int64_t sum_dim_up = 0;
      for (auto [name, mult] : rep_name_mult) {
        auto irrep = read_representation(fl, name);
        auto block = Spinhalf(nsites, nup, irrep);
        sum_dim_up += block.size() * mult;
        sum_dim += block.size() * mult;
      }
      REQUIRE(sum_dim_up == math::binomial(nsites, nup));
    }
    REQUIRE(sum_dim == ((int64_t)1 << nsites));
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: apply correctness in symmetry-resolved sectors
// -----------------------------------------------------------------------

TEST_CASE("spinhalf_symmetric_apply", "[spinhalf]") try {

  // Cyclic Heisenberg chains with J1=J2=1
  for (int64_t nsites = 3; nsites < 7; ++nsites) {
    Log("spinhalf_symmetric_apply: HB chain, N={}", nsites);
    std::vector<Representation> irreps;
    for (int64_t k = 0; k < nsites; ++k) {
      irreps.push_back(cyclic_group_irrep(nsites, k));
    }
    auto ops = HBchain(nsites, 1.0, 1.0);
    test_spinhalf_symmetric_apply(ops, nsites, irreps);
    test_spinhalf_symmetric_apply_no_sz(ops, nsites, irreps);
  }

  // Triangular 3x3
  {
    Log("spinhalf_symmetric_apply: triangular 3x3");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["Jz1"] = 1.00;
    ops["Jz2"] = 0.23;
    ops["Jx1"] = 0.76;
    ops["Jx2"] = 0.46;

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
        {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
        {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
        {"Y.D1.A", 6},      {"Y.D1.B", 6}};

    std::vector<Representation> irreps;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(fl, name));
    }
    test_spinhalf_symmetric_apply(ops, 9, irreps);
    test_spinhalf_symmetric_apply_no_sz(ops, 9, irreps);
  }

  // Triangular J1-J2-Jchi N=12: check ground state energy per irrep
  {
    Log("spinhalf_symmetric_apply: triangular J1J2Jchi N=12 ground states");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = -0.09;

    std::vector<std::pair<std::string, double>> rep_name_e0 = {
        {"Gamma.C6.A", -6.945600070082439181},
        {"Gamma.C6.B", -5.84109124378730904},
	// used wrong character conjugation
        // {"Gamma.C6.E1a", -3.845414358083013795},
        // {"Gamma.C6.E1b", -3.855641724835592754},
        // {"Gamma.C6.E2a", -6.322314852895261517},
        // {"Gamma.C6.E2b", -6.415724335805900402},
        {"Gamma.C6.E1a", -3.855641724835592754},
        {"Gamma.C6.E1b", -3.845414358083013795},
        {"Gamma.C6.E2a", -6.415724335805900402},
        {"Gamma.C6.E2b", -6.322314852895261517},
        {"K.C3.A", -5.91975118116224408},
        // {"K.C3.Ea", -5.204513347364086329},
        // {"K.C3.Eb", -5.028170383600089721},
	{"K.C3.Ea", -5.028170383600089721},
        {"K.C3.Eb", -5.204513347364086329},
        {"M.C2.A", -5.756684675081961799},
        {"M.C2.B", -5.77235103255616977},
        {"X.C1.A", -5.90306276605228053}};

    int64_t nsites = 12;
    int64_t nup = 6;
    for (auto [name, energy] : rep_name_e0) {
      auto irrep = read_representation(fl, name);
      auto spinhalf = Spinhalf(nsites, nup, irrep);
      double e0 = eigval0(ops, spinhalf);
      REQUIRE(isapprox(e0, energy, 1e-12, 1e-10));
    }
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: symmetry-resolved spectra reconstruct the full spectrum
// -----------------------------------------------------------------------

TEST_CASE("spinhalf_symmetric_matrix", "[spinhalf]") try {

  // Cyclic Heisenberg chains with J1=J2=1
  for (int64_t nsites = 3; nsites < 7; ++nsites) {
    Log("spinhalf_symmetric_matrix: HB chain, N={}", nsites);
    std::vector<Representation> irreps;
    std::vector<int64_t> multiplicities;
    for (int64_t k = 0; k < nsites; ++k) {
      irreps.push_back(cyclic_group_irrep(nsites, k));
      multiplicities.push_back(1);
    }
    auto ops = HBchain(nsites, 1.0, 1.0);
    test_spinhalf_symmetric_spectra(ops, nsites, irreps, multiplicities);
    test_spinhalf_symmetric_spectra_no_sz(ops, nsites, irreps, multiplicities);
  }

  // Triangular 3x3
  {
    Log("spinhalf_symmetric_matrix: triangular 3x3");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["Jz1"] = 1.00;
    ops["Jz2"] = 0.23;
    ops["Jx1"] = 0.76;
    ops["Jx2"] = 0.46;

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
        {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
        {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
        {"Y.D1.A", 6},      {"Y.D1.B", 6}};

    std::vector<Representation> irreps;
    std::vector<int64_t> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(fl, name));
      multiplicities.push_back(mult);
    }
    test_spinhalf_symmetric_spectra(ops, 9, irreps, multiplicities);
    test_spinhalf_symmetric_spectra_no_sz(ops, 9, irreps, multiplicities);
  }

  // Triangular J1-J2-Jchi N=12: verify ground state energy per irrep via full
  // diag
  {
    Log("spinhalf_symmetric_matrix: triangular J1J2Jchi N=12");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = -0.09;

    std::vector<std::pair<std::string, double>> rep_name_e0 = {
        {"Gamma.C6.A", -6.945600070082439181},
        {"Gamma.C6.B", -5.84109124378730904},
	// used wrong character conjugation
        // {"Gamma.C6.E1a", -3.845414358083013795},
        // {"Gamma.C6.E1b", -3.855641724835592754},
        // {"Gamma.C6.E2a", -6.322314852895261517},
        // {"Gamma.C6.E2b", -6.415724335805900402},
        {"Gamma.C6.E1a", -3.855641724835592754},
        {"Gamma.C6.E1b", -3.845414358083013795},
        {"Gamma.C6.E2a", -6.415724335805900402},
        {"Gamma.C6.E2b", -6.322314852895261517},
        {"K.C3.A", -5.91975118116224408},
        // {"K.C3.Ea", -5.204513347364086329},
        // {"K.C3.Eb", -5.028170383600089721},
	{"K.C3.Ea", -5.028170383600089721},
        {"K.C3.Eb", -5.204513347364086329},
        {"M.C2.A", -5.756684675081961799},
        {"M.C2.B", -5.77235103255616977},
        {"X.C1.A", -5.90306276605228053}};

    int64_t nsites = 12;
    int64_t nup = 6;
    for (auto [name, energy] : rep_name_e0) {
      auto irrep = read_representation(fl, name);
      auto spinhalf = Spinhalf(nsites, nup, irrep);
      auto H = matrixC(ops, spinhalf, spinhalf);
      REQUIRE(arma::norm(H - H.t()) < 1e-12);
      arma::vec eigs;
      arma::eig_sym(eigs, H);
      REQUIRE(isapprox(eigs(0), energy, 1e-12, 1e-10));
    }
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

// -----------------------------------------------------------------------
// TEST_CASE: Kitaev-Gamma model on N=8 honeycomb
// -----------------------------------------------------------------------

TEST_CASE("kitaev_gamma", "[spinhalf]") {

  Log("kitaev_gamma: Kitaev-Gamma honeycomb N=8");

  {
    double K = -.30901699437494742412;
    double G = .95105651629515357210;
    test_kitaev_gamma(K, G,
                      {{"Gamma.C2.A", -3.1765766652975568896},
                       {"Gamma.C2.B", -2.0558245985778302867},
                       {"M0.C2.A", -2.3210367200223522843},
                       {"M0.C2.B", -2.8044982600737728973},
                       {"M1.C2.A", -2.3210367200223518402},
                       {"M1.C2.B", -2.8044982600737733414},
                       {"M2.C2.A", -2.3210367200223536166},
                       {"M2.C2.B", -2.8044982600737742295}});
  }

  {
    double K = -1.0;
    double G = 0.0;
    test_kitaev_gamma(K, G,
                      {{"Gamma.C2.A", -1.732050807568876083},
                       {"Gamma.C2.B", -0.85463767971846216209},
                       {"M0.C2.A", -1.4516059629557762634},
                       {"M0.C2.B", -1.4516059629557771515},
                       {"M1.C2.A", -1.4516059629557773736},
                       {"M1.C2.B", -1.4516059629557787058},
                       {"M2.C2.A", -1.4516059629557762634},
                       {"M2.C2.B", -1.4516059629557764854}});
  }

  {
    double K = 0.0;
    double G = 1.0;
    test_kitaev_gamma(K, G,
                      {{"Gamma.C2.A", -2.9261296954884912225},
                       {"Gamma.C2.B", -2.2611516188994311705},
                       {"M0.C2.A", -2.3314966606292610862},
                       {"M0.C2.B", -2.6582065612437091318},
                       {"M1.C2.A", -2.3314966606292633067},
                       {"M1.C2.B", -2.6582065612437086877},
                       {"M2.C2.A", -2.3314966606292646389},
                       {"M2.C2.B", -2.6582065612437109081}});
  }

  {
    double K = -0.64944804833018365819;
    double G = 0.76040596560003093085;
    test_kitaev_gamma(K, G,
                      {{"Gamma.C2.A", -3.0857997869548752234},
                       {"Gamma.C2.B", -1.8283034696482216575},
                       {"M0.C2.A", -2.0356565020285248835},
                       {"M0.C2.B", -2.6544975593746693576},
                       {"M1.C2.A", -2.035656502028527548},
                       {"M1.C2.B", -2.6544975593746698017},
                       {"M2.C2.A", -2.0356565020285266598},
                       {"M2.C2.B", -2.6544975593746689135}});
  }

  {
    double K = -0.70710678118654752441;
    double G = 0.70710678118654752441;
    test_kitaev_gamma(K, G,
                      {{"Gamma.C2.A", -3.0140686352095316103},
                       {"Gamma.C2.B", -1.7541857925174297872},
                       {"M0.C2.A", -1.94808553693512998},
                       {"M0.C2.B", -2.5821356185598034472},
                       {"M1.C2.A", -1.9480855369351295359},
                       {"M1.C2.B", -2.5821356185598021149},
                       {"M2.C2.A", -1.948085536935130424},
                       {"M2.C2.B", -2.582135618559802559}});
  }

  Log("kitaev_gamma: done");
}

// -----------------------------------------------------------------------
// TEST_CASE: test whether exchange terms constructed correctly
// -----------------------------------------------------------------------

TEST_CASE("spinhalfexchange", "[spinhalf]") try {

  auto b = Spinhalf(2);
  {
    auto op1 = complex(0.0, 1.0) * Op("Exchange", {0, 1});
    auto op2 = complex(0.0, 1.0) * Op("Exchange", {1, 0});

    // XDIAG_SHOW(op1);
    // XDIAG_SHOW(hc(op1));

    // auto algebra = algebra::matrix_algebra(2);
    auto algebra = algebra::spin_algebra(2);

    // XDIAG_SHOW(normal_order(op1, algebra));
    // XDIAG_SHOW(normal_order(op2, algebra));
    REQUIRE(isapprox(op1, op2, algebra));
    REQUIRE_FALSE(isapprox(hc(op1), op1, algebra));
    REQUIRE_FALSE(isapprox(hc(op2), op1, algebra));

    auto m1 = matrixC(op1, b);
    auto m2 = matrixC(op2, b);
    REQUIRE(isapprox(m1, m2));
  }
  {
    auto op1 = complex(0.0, 1.0) * Op("ExchangeAsym", {0, 1});
    auto op2 = complex(0.0, 1.0) * Op("ExchangeAsym", {1, 0});

    // XDIAG_SHOW(op1);
    // XDIAG_SHOW(hc(op1));

    // auto algebra = algebra::matrix_algebra(2);
    auto algebra = algebra::spin_algebra(2);

    // XDIAG_SHOW(normal_order(op1, algebra));
    // XDIAG_SHOW(normal_order(op2, algebra));
    REQUIRE(isapprox(op1, -op2, algebra));
    REQUIRE(isapprox(hc(op1), op1, algebra));
    REQUIRE(isapprox(hc(op2), op2, algebra));

    auto m1 = matrixC(op1, b);
    auto m2 = matrixC(op2, b);
    REQUIRE(isapprox(m1, arma::cx_mat(-m2)));
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

// Randomized cross-check of the full operator pipeline against naive matrix
// products (shared harness). Uses the full spin-1/2 Hilbert space (no Sz
// conservation) so every elementary op is an endomorphism.
TEST_CASE("spinhalfrandomopsum", "[spinhalf]") try {
  for (int nsites = 3; nsites < 6; ++nsites) {
    Log("Spinhalf random OpSum matrix test: N = {}", nsites);
    for (uint32_t seed = 0; seed < 5; ++seed) {
      testcases::test_random_opsum_matrix(Spinhalf(nsites), seed);
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
