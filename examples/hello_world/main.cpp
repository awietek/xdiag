#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  say_hello();

  auto block = Electron(2, 1, 1);
  auto bond = Bond("HB", 1.0, {0});
  auto mat = matrix(bond, block);
  XDIAG_PRINT(mat);
} catch (std::exception const &e) {
  traceback(e);
}
