#include "../catch.hpp"

#include <extern/armadillo/armadillo>

#include <hydra/common.h>
#include <hydra/io/args.h>

TEST_CASE("args", "[io]") {
  using namespace hydra;
  using namespace arma;

  Args args({{"bool", true},
             {"string", "hello"},
             {"int", 42},
             {"double", 1.23},
             {"complex", 1.23 + 3.21i}});

  REQUIRE(args.defined("bool"));
  REQUIRE(args.defined("string"));
  REQUIRE(args.defined("int"));
  REQUIRE(args.defined("double"));
  REQUIRE(args.defined("complex"));
  REQUIRE(!args.defined("vec"));

  REQUIRE(args["bool"].as<bool>() == true);
  REQUIRE(args["string"].as<std::string>() == "hello");
  REQUIRE(args["int"].as<int>() == 42);
  REQUIRE(args["double"].as<double>() == 1.23);
  REQUIRE(args["complex"].as<complex>() == 1.23 + 3.21i);

  REQUIRE(args["asdf0"].as<bool>(false) == false);
  REQUIRE(args["asdf1"].as<std::string>("jjj") == "jjj");
  REQUIRE(args["asdf2"].as<int>(43) == 43);
  REQUIRE(args["asdf3"].as<double>(4.3) == 4.3);
  REQUIRE(args["asdf4"].as<complex>(4.3 + 1.2i) == 4.3 + 1.2i);
}
