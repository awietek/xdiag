#include "../catch.hpp"

#include <hydra/all.h>

TEST_CASE("file_toml", "[io]") {
  using namespace hydra;

  // Just try to parse everything in the example toml file
  std::string filename = "data/toml/test.toml";
  auto fl = FileToml(filename, 'r');

  // Parse String
  std::string n = "title";
  Log("{}", fl[n].as<std::string>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<std::string>() == "TOML Beispiel");

  // Parse integers
  n = "datenbank.verbindungen_max";
  Log("{}", fl[n].as<int32_t>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<int8_t>() == 42);
  REQUIRE(fl[n].as<int16_t>() == 42);
  REQUIRE(fl[n].as<int32_t>() == 42);
  REQUIRE(fl[n].as<int64_t>() == 42);
  REQUIRE(fl[n].as<uint8_t>() == 42);
  REQUIRE(fl[n].as<uint16_t>() == 42);
  REQUIRE(fl[n].as<uint32_t>() == 42);
  REQUIRE(fl[n].as<uint64_t>() == 42);

  // Parse bool
  n = "datenbank.aktiviert";
  Log("{}", fl[n].as<bool>());
  std::cout << fl[n].as<bool>() << "\n";
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<bool>() == true);

  n = "floating.pi";
  Log("{}", fl[n].as<double>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<double>() == 3.1415);

  n = "floating.cplx";
  Log("{}", fl[n].as<complex>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<complex>() == 2.3122 + .1237i);
    
}
