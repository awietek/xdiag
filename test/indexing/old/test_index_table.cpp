#include "../catch.hpp"

#include <hydra/all.h>

template <class basis_t, class idx_t=uint32>
void test_index_table(basis_t const &basis) {
  using namespace hydra;

  IndexTable<basis_t> index_table(basis);

  idx_t idx = 0;
  for (auto state : basis) {
    // std::cout << String(state, basis.n_sites()) << "\n";
    REQUIRE(index_table.index(state) == idx);
    REQUIRE(index_table.state(idx) == state);
    ++idx;
  }
}

TEST_CASE("indexing/index_table", "[indexing]") {
  using namespace hydra;

  // Spinhalf
  for (int n_sites = 0; n_sites < 7; ++n_sites)
    for (int n_up = 0; n_up <= n_sites; ++n_up) {
      BasisSpinHalf<uint32> basis(n_sites, n_up);
      test_index_table(basis);
    }

  // TJ
  for (int n_sites = 0; n_sites < 6; ++n_sites)
    for (int n_up = 0; n_up <= n_sites; ++n_up)
      for (int n_dn = 0; n_dn <= n_sites - n_up; ++n_dn) {
        BasisTJ<uint32> basis(n_sites, {n_up, n_dn});
        test_index_table(basis);
      }

  // Electron
  for (int n_sites = 0; n_sites < 5; ++n_sites)
    for (int n_up = 0; n_up <= n_sites; ++n_up)
      for (int n_dn = 0; n_dn <= n_sites; ++n_dn) {
        BasisElectron<uint32> basis(n_sites, {n_up, n_dn});
  	test_index_table(basis);
      }	
}
