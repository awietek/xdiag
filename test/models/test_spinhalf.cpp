#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;

template <class bit_t>
void print_states(Spinhalf<bit_t> spinhalf){
  int n_sites = spinhalf.n_sites();
  for (bit_t state : spinhalf){
    std::cout << bits_to_string(state, n_sites) << "\n";
  }
}

template <class bit_t> void test_spinhalf(int n_sites) {
  
  // Full spinhalf space
  std::cout << "Sz: nA\n";
  auto spinhalf1 = Spinhalf<bit_t>(n_sites);
  print_states(spinhalf1);
  std::cout <<  "\n";

  // spinhalf space sz conserved
  for (int sz=-n_sites/2; sz<=n_sites/2; ++sz){
    std::cout << "Sz: " << sz << "\n";
    auto spinhalf2 = Spinhalf<bit_t>(n_sites, sz);
    print_states(spinhalf2);
    std::cout <<  "\n";
  }

  // cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = SpaceGroup<bit_t>(permutations);
  auto irrep = Representation({"Gamma", std::vector<complex>(n_sites, 1.0)});

  // Space group symmetries / no Sz conservation
  std::cout << "Sz: nA cyclic\n";
  auto spinhalf3 = Spinhalf<bit_t>(n_sites, space_group, irrep);
  print_states(spinhalf3);
  std::cout <<  "\n";

  // Space group symmetries / Sz conserved
  for (int sz=-n_sites/2; sz<=n_sites/2; ++sz){
    std::cout << "Sz: " << sz << " cyclic\n";
    auto spinhalf4 = Spinhalf<bit_t>(n_sites, sz, space_group, irrep);
    print_states(spinhalf4);
    std::cout <<  "\n";
  }

}

TEST_CASE("spinhalf", "[models]") {
  test_spinhalf<uint32>(4);
  
  // for (int n_sites=1; n_sites<6; ++n_sites) {
  //   test_spinhalf<uint16>(n_sites);
  //   test_spinhalf<uint32>(n_sites);
  //   test_spinhalf<uint64>(n_sites);
  // }

}
