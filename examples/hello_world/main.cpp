#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  say_hello();
} catch (Error e) {
  error_trace(e);
}
