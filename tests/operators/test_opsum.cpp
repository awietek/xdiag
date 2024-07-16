#include "../catch.hpp"

#include <iostream>
#include <xdiag/operators/opsum.hpp>

TEST_CASE("opsum", "[operators]") {
  using namespace xdiag;

  OpSum ops;
  ops += Op("HB", "J1", {0, 1});
  ops += Op("HB", "J1", {1, 2});
  ops += Op("HB", "J1", {2, 3});
  ops += Op("HB", "J1", {3, 4});
  ops += Op("HB", "J1", {4, 5});
  ops += Op("HB", "J1", {5, 0});
  ops += Op("HB", "J2", {0, 2});
  ops += Op("HB", "J2", {1, 3});
  ops += Op("HB", "J2", {2, 4});
  ops += Op("HB", "J2", {3, 5});
  ops += Op("HB", "J2", {4, 0});
  ops += Op("HB", "J2", {5, 1});
  {
    Op b1("A", "B", {1, 2, 3});
    Op b2("A", "B", {2, 3, 4});
    std::vector<int64_t> c = {2, 3};
    REQUIRE(common_sites(b1, b2) == c);
  }
  {
    Op b1("A", "B", {1, 2, 3});
    Op b2("A", "B", {6, 7, 8});
    std::vector<int64_t> c = {};
    REQUIRE(common_sites(b1, b2) == c);
  }
  {
    Op b1("A", "B", {2, 2, 3});
    Op b2("A", "B", {2, 7, 8});
    std::vector<int64_t> c = {2};
    REQUIRE(common_sites(b1, b2) == c);
  }

  // for (auto op : ops)
  //   std::cout += op += std::endl;
  // std::cout += std::endl;
  // for (auto op : ops.ops_of_coupling("J1"))
  //   std::cout += op += std::endl;
  // std::cout += std::endl;
  // for (auto op : ops.ops_of_coupling("J2"))
  //   std::cout += op += std::endl;
}
