#include "../catch.hpp"

#include "../blocks/electron/testcases_electron.h"

#include <xdiag/symmetries/generated_group.h>
#include <xdiag/utils/close.h>

TEST_CASE("generated_group", "[symmetries]") {
  using namespace xdiag;
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;

  Log("Testing generated_group for linear chain");
  for (int64_t n_sites = 3; n_sites < 7; ++n_sites) {
    auto [group, irreps, multiplicities] =
        get_cyclic_group_irreps_mult(n_sites);
    (void)multiplicities;

    std::vector<int64_t> translation;
    for (int64_t s = 0; s < n_sites; ++s) {
      translation.push_back((s + 1) % n_sites);
    }
    Permutation perm(translation);

    auto group2 = generated_group(perm);
    REQUIRE(group == group2);

    for (int64_t k = 0; k < n_sites; ++k) {
      auto irrep = irreps[k];

      std::vector<complex> chis;
      for (int64_t l = 0; l < n_sites; ++l) {
        chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                        std::sin(2 * M_PI * l * k / n_sites)});
      }
      auto irrep2 = Representation(chis);
      REQUIRE(irrep == irrep2);
    }
  }
  Log("done");

  {
    Log("Testing generated_group for 2x2 square");
    Permutation t1{1, 0, 3, 2};
    Permutation t2{2, 3, 0, 1};
    auto group = generated_group({t1, t2});
    REQUIRE(group.size() == 4);

    auto irrep1 = generated_irrep({t1, t2}, {1.0, 1.0});
    REQUIRE(irrep1.characters() == std::vector<complex>{1.0, 1.0, 1.0, 1.0});

    auto irrep2 = generated_irrep({t1, t2}, {1.0, -1.0});
    REQUIRE(irrep2.characters() == std::vector<complex>{1.0, 1.0, -1.0, -1.0});

    auto irrep3 = generated_irrep({t1, t2}, {-1.0, 1.0});
    REQUIRE(irrep3.characters() == std::vector<complex>{1.0, -1.0, 1.0, -1.0});

    auto irrep4 = generated_irrep({t1, t2}, {-1.0, -1.0});
    REQUIRE(irrep4.characters() == std::vector<complex>{1.0, -1.0, -1.0, 1.0});
    Log("done");
  }

  {
    Log("Testing generated_group for 3x3 square");
    Permutation t1{1, 2, 0, 4, 5, 3, 7, 8, 6};
    Permutation t2{3, 4, 5, 6, 7, 8, 0, 1, 2};
    auto group = generated_group({t1, t2});
    REQUIRE(group.size() == 9);

    std::vector<complex> chis = {
        1.0,
        {std::cos(2 * M_PI / 3), std::sin(2 * M_PI / 3)},
        {std::cos(4 * M_PI / 3), std::sin(4 * M_PI / 3)}};

    for (auto c1 : chis) {
      for (auto c2 : chis) {
        auto irrep = generated_irrep({t1, t2}, {c1, c2});
        arma::cx_vec chars = {1.0,          c1,           c2,
                              c1 * c1,      c1 * c2,      c2 * c2,
                              c1 * c1 * c2, c1 * c2 * c2, c1 * c1 * c2 * c2};
        REQUIRE(close(arma::cx_vec(irrep.characters()), chars));
      }
    }
    Log("done");
  }
}
