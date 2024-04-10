#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  say_hello();
} catch (std::exception const &e) {
  traceback(e);
}
