#include "../../catch.hpp"

#include <iostream>

using namespace xdiag;

template <typename bit_t>
void test_electron_symmetric_indexing(PermutationGroup const &group,
                                      Representation const &irrep) {
  int n_sites = group.n_sites();
  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      auto idxer = ElectronSymmetricIndexing(n_sites, nup, ndn, group, irrep)
    }
  }
}

TEST_CASE("electron_symmetric_indexing", "[indexing]") {
  Log("Test ElectronSymmetricIndexing");
  int max_N = 6;

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    Log("chain N={}", n_sites);
    auto [group, irreps] =
        xdiag::testcases::electron::get_cyclic_group_irreps(n_sites);
    test_fermi_bool_table<uint16_t>(group);
    test_fermi_bool_table<uint32_t>(group);
    test_fermi_bool_table<uint64_t>(group);
  }

  Log("3x3 triangular");
  using bit_t = uint16_t;
  std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";
  auto permutations = xdiag::read_permutations(lfile);
  auto group = PermutationGroup(permutations);

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};
  irreps.clear();
  multiplicities.clear();
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_represenation(lfile, name));
    multiplicities.push_back(mult);
  }
  test_electron_symmetric_spectra<bit_t>(ops, couplings, space_group,
                                         irreps, multiplicities);
  Log("done");
}
