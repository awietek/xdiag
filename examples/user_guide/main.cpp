
#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
// clang-format off
 
{
// --8<-- [start:first_steps_1]
using namespace xdiag;
int N = 8;
auto hspace = Spinhalf(N);
// --8<-- [end:first_steps_1]

// --8<-- [start:first_steps_2]
for (auto spins : hspace) {
  Log("{}", to_string(spins));
}
// --8<-- [end:first_steps_2]

// --8<-- [start:first_steps_3]
int nup = 4;
auto block = Spinhalf(N, nup);
for (auto spins : block) {
  Log("{}", to_string(spins));
}
// --8<-- [end:first_steps_3]

// --8<-- [start:first_steps_4]
XDIAG_SHOW(hspace.size());
XDIAG_SHOW(block.size());
// --8<-- [end:first_steps_4]

 
}
  
// clang-format on

} catch (Error e) {
  error_trace(e);
}
