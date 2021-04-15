#include "../catch.hpp"

#include <hydra/all.h>

#include <chrono>
#include <iostream>

TEST_CASE( "binomial", "[combinatorics]" ) {
  using namespace hydra;
  using namespace hydra::combinatorics;

  // Check Pascals triangle manually
  REQUIRE(binomial(1, -1) == 0);
  REQUIRE(binomial(1, 0) == 1);
  REQUIRE(binomial(1, 1) == 1);
  REQUIRE(binomial(1, 2) == 0);

  REQUIRE(binomial(2, -1) == 0);
  REQUIRE(binomial(2, 0) == 1);
  REQUIRE(binomial(2, 1) == 2);
  REQUIRE(binomial(2, 2) == 1);
  REQUIRE(binomial(2, 3) == 0);

  REQUIRE(binomial(3, -1) == 0);
  REQUIRE(binomial(3, 0) == 1);
  REQUIRE(binomial(3, 1) == 3);
  REQUIRE(binomial(3, 2) == 3);
  REQUIRE(binomial(3, 3) == 1);
  REQUIRE(binomial(3, 4) == 0);

  REQUIRE(binomial(4, -1) == 0);
  REQUIRE(binomial(4, 0) == 1);
  REQUIRE(binomial(4, 1) == 4);
  REQUIRE(binomial(4, 2) == 6);
  REQUIRE(binomial(4, 3) == 4);
  REQUIRE(binomial(4, 4) == 1);
  REQUIRE(binomial(4, 5) == 0);

  REQUIRE(binomial(5, -1) == 0);
  REQUIRE(binomial(5, 0) == 1);
  REQUIRE(binomial(5, 1) == 5);
  REQUIRE(binomial(5, 2) == 10);
  REQUIRE(binomial(5, 3) == 10);
  REQUIRE(binomial(5, 4) == 5);
  REQUIRE(binomial(5, 5) == 1);
  REQUIRE(binomial(5, 6) == 0);

  REQUIRE(binomial(6, -1) == 0);
  REQUIRE(binomial(6, 0) == 1);
  REQUIRE(binomial(6, 1) == 6);
  REQUIRE(binomial(6, 2) == 15);
  REQUIRE(binomial(6, 3) == 20);
  REQUIRE(binomial(6, 4) == 15);
  REQUIRE(binomial(6, 5) == 6);
  REQUIRE(binomial(6, 6) == 1);
  REQUIRE(binomial(6, 7) == 0);

  // Check whether lookups give sam as computation
  for (int n=0; n<40; ++n)
    for (int k=0; k<40; ++k)
      REQUIRE(binomial(n, k) == binom(n, k));

  // // time lookup vs computation of binomials
  // using clock = std::chrono::high_resolution_clock;
  // using secs = std::chrono::duration<double>;

  // int64 sum = 0;
  // auto t1 = clock::now();
  // for (int i=0; i<10000; ++i)
  //   for (int n=0; n<32; ++n)
  //     for (int k=0; k<32; ++k)
  // 	{
  // 	  sum += binomial(n, k);
  // 	}
  // auto t2 = clock::now();
  // std::cout << "lookup " << secs(t2-t1).count() << " " << sum << "\n";

  // sum = 0;
  // t1 = clock::now();
  // for (int i=0; i<10000; ++i)
  //   for (int n=0; n<32; ++n)
  //     for (int k=0; k<32; ++k)
  // 	{
  // 	  sum += binom(n, k);
  // 	}
  // t2 = clock::now();
  // std::cout << "compute " << secs(t2-t1).count() << " " << sum << "\n";

}
