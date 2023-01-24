#include "../catch.hpp"

#include <hydra/all.h>

TEST_CASE("file_toml", "[io]") {
  using namespace hydra;

  std::string filename = "data/toml/test.toml";
  auto fl = FileToml(filename, 'r');

  REQUIRE(fl.defined("title") &&
          (fl["title"].as<std::string>() == "TOML Beispiel"));
}
