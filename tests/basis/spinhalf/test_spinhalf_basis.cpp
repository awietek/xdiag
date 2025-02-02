#include "../../catch.hpp"

#include <iostream>

#include <iostream>
#include <xdiag/basis/spinhalf/basis_symmetric_no_sz.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_sz.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/vector.hpp>

using namespace xdiag;

template <typename bit_t, class Basis>
void test_spinhalf_basis_state(Basis const &basis, bit_t state) {

  auto const &group_action = basis.group_action();
  auto const &irrep = basis.irrep();

  int64_t idx_rep = basis.index(state);
  if (idx_rep == invalid_index) {
    if (isreal(irrep)) {
      Vector chs = irrep.characters();
      double norm = symmetries::norm(state, group_action, chs.as<arma::vec>());
      REQUIRE(std::abs(norm) < 1e-6);
    } else {
      Vector chs = irrep.characters();
      double norm =
          symmetries::norm(state, group_action, chs.as<arma::cx_vec>());
      REQUIRE(std::abs(norm) < 1e-6);
    }
  } else {

    bit_t rep = symmetries::representative(state, group_action);
    bit_t rep_lookup = basis.representative(state);
    REQUIRE(rep == rep_lookup);

    REQUIRE(rep == basis.state(idx_rep));

    // int nsites = basis.nsites();
    // std::cout << BSTR(state) << "\n";

    // double norm = symmetries::norm(state, group_action, irrep);
    // REQUIRE(norm == basis.norm(idx_rep));

    if (isreal(irrep)) {
      Vector chs = irrep.characters();
      double norm = symmetries::norm(state, group_action, chs.as<arma::vec>());
      REQUIRE(norm == basis.norm(idx_rep));
    } else {
      Vector chs = irrep.characters();
      double norm =
          symmetries::norm(state, group_action, chs.as<arma::cx_vec>());
      REQUIRE(norm == basis.norm(idx_rep));
    }

    auto [idxr1, sym] = basis.index_sym(state);
    (void)idxr1;
    REQUIRE(group_action.apply(sym, state) == rep);

    auto [idxr2, syms] = basis.index_syms(state);
    (void)idxr2;
    for (int sym : syms) {
      REQUIRE(group_action.apply(sym, state) == rep);
    }
  }
}

template <typename bit_t>
void check_basis_symmetric_sz(
    basis::spinhalf::BasisSymmetricSz<bit_t> const &basis) {

  int nsites = basis.nsites();
  int nup = basis.nup();

  if (basis.size() > 0) {
    for (auto state : combinatorics::Combinations<uint64_t>(nsites, nup)) {
      test_spinhalf_basis_state(basis, state);
    }
  }
}

template <typename bit_t>
void check_basis_symmetric_no_sz(
    basis::spinhalf::BasisSymmetricNoSz<bit_t> const &basis) {

  int nsites = basis.nsites();
  if (basis.size() > 0) {

    for (auto state : combinatorics::Subsets<uint64_t>(nsites)) {
      test_spinhalf_basis_state(basis, state);
    }
  }
}

template <class bit_t> void test_basis_symmetric() {

  using namespace xdiag::basis::spinhalf;

  Log("BasisSymmetric: triangular 3x3");
  int nsites = 9;
  std::string lfile = XDIAG_DIRECTORY
      "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
  auto fl = FileToml(lfile);
  // std::vector<std::string> irrep_names = {
  //     "Gamma.D6.A1", "Gamma.D6.A2", "Gamma.D6.B1", "Gamma.D6.B2",
  //     "Gamma.D6.E1", "Gamma.D6.E2", "K.D3.A1",     "K.D3.A2",
  //     "K.D3.E",      "Y.D1.A",      "Y.D1.B"};

  std::vector<std::string> irrep_names = {"Gamma.D6.B2"};

  for (auto irrep_name : irrep_names) {
    Log("irrep: {}", irrep_name);
    auto irrep = read_representation(fl, irrep_name);
    auto idxng = BasisSymmetricNoSz<bit_t>(irrep);
    // check_basis_symmetric_no_sz<bit_t>(idxng);

    for (int nup = 3; nup <= 3; ++nup) {
      auto idxng = BasisSymmetricSz<bit_t>(nup, irrep);
      check_basis_symmetric_sz<bit_t>(idxng);
    }
  }
}

TEST_CASE("spinhalf_basis", "[basis]") {
  Log("Test SpinhalfBasisSublattice");
  Log("uint32_t");
  test_basis_symmetric<uint32_t>();
  Log("uint64_t");
  test_basis_symmetric<uint64_t>();
  Log("Done");
}
