#include <hydra/all.h>

using namespace hydra;

int main() try {
  say_hello();
} catch (std::exception const &e) {
  traceback(e);
}
