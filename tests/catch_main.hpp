#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include <xdiag/utils/error.hpp>

int main(int argc, char *argv[]) try {
  int result = Catch::Session().run(argc, argv);
  return result;
} catch (xdiag::Error const &e) {
  error_trace(e);
}
