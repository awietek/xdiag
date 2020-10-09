#include "catch.hpp"

#include <hydra/all.h>

TEST_CASE("qns", "[qns]") {
  using namespace hydra;

  int min_n_sites = 4;
  int max_n_sites = 6;

  // Test qn_spin
  for (int n_sites = min_n_sites; n_sites < max_n_sites; ++n_sites) {

    for (int n_up1 = -n_sites; n_up1 < 2 * n_sites; ++n_up1) {
      qn_spinhalf q1 = {n_up1};
      if ((0 <= n_up1) && (n_up1 <= n_sites))
        REQUIRE(valid(q1, n_sites));
      else
        REQUIRE(!valid(q1, n_sites));

      for (int n_up2 = -n_sites; n_up2 < 2 * n_sites; ++n_up2) {
        qn_spinhalf q2 = {n_up2};

        auto q_add = q1 + q2;
        REQUIRE(q_add.n_up == q1.n_up + q2.n_up);

        auto q_sub = q1 - q2;
        REQUIRE(q_sub.n_up == q1.n_up - q2.n_up);

        if (n_up1 == n_up2)
          REQUIRE(q1 == q2);
        else
          REQUIRE(q1 != q2);
      }
    }
  }

  // Test qn_electron
  for (int n_sites = min_n_sites; n_sites < max_n_sites; ++n_sites) {

    for (int n_up1 = -n_sites; n_up1 < 2 * n_sites; ++n_up1)
      for (int n_dn1 = -n_sites; n_dn1 < 2 * n_sites; ++n_dn1) {

        qn_electron q1 = {n_up1, n_dn1};
        if (((0 <= n_up1) && (n_up1 <= n_sites)) &&
            ((0 <= n_dn1) && (n_dn1 <= n_sites))) {
          REQUIRE(valid(q1, n_sites));
        } else
          REQUIRE(!valid(q1, n_sites));

        for (int n_up2 = -n_sites; n_up2 < 2 * n_sites; ++n_up2)
          for (int n_dn2 = -n_sites; n_dn2 < 2 * n_sites; ++n_dn2) {
            {
              qn_electron q2 = {n_up2, n_dn2};

              auto q_add = q1 + q2;
              REQUIRE(q_add.n_up == q1.n_up + q2.n_up);
              REQUIRE(q_add.n_dn == q1.n_dn + q2.n_dn);

              auto q_sub = q1 - q2;
              REQUIRE(q_sub.n_up == q1.n_up - q2.n_up);
              REQUIRE(q_sub.n_dn == q1.n_dn - q2.n_dn);

              if ((n_up1 == n_up2) && (n_dn1 == n_dn2))
                REQUIRE(q1 == q2);
              else
                REQUIRE(q1 != q2);
            }
          }
      }
  }

  // Test qn_tj
  for (int n_sites = min_n_sites; n_sites < max_n_sites; ++n_sites) {

    for (int n_up1 = -n_sites; n_up1 < 2 * n_sites; ++n_up1)
      for (int n_dn1 = -n_sites; n_dn1 < 2 * n_sites; ++n_dn1) {

        qn_tj q1 = {n_up1, n_dn1};
        if (((0 <= n_up1) && (n_up1 <= n_sites)) &&
            ((0 <= n_dn1) && (n_dn1 <= n_sites)) &&
            (n_up1 + n_dn1 <= n_sites)) {
          REQUIRE(valid(q1, n_sites));
        } else
          REQUIRE(!valid(q1, n_sites));

        for (int n_up2 = -n_sites; n_up2 < 2 * n_sites; ++n_up2)
          for (int n_dn2 = -n_sites; n_dn2 < 2 * n_sites; ++n_dn2) {
            {
              qn_tj q2 = {n_up2, n_dn2};

              auto q_add = q1 + q2;
              REQUIRE(q_add.n_up == q1.n_up + q2.n_up);
              REQUIRE(q_add.n_dn == q1.n_dn + q2.n_dn);

              auto q_sub = q1 - q2;
              REQUIRE(q_sub.n_up == q1.n_up - q2.n_up);
              REQUIRE(q_sub.n_dn == q1.n_dn - q2.n_dn);

              if ((n_up1 == n_up2) && (n_dn1 == n_dn2))
                REQUIRE(q1 == q2);
              else
                REQUIRE(q1 != q2);
            }
          }
      }
  }
}
