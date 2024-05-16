#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  // clang-format off

{
// --8<-- [start:Permutation]
Permutation p1 = {0, 2, 1, 3};
Permutation p2 = {2, 0, 1, 3};

XDIAG_PRINT(inverse(p1));
XDIAG_PRINT(p1*p2);
// --8<-- [end:Permutation]
}

{
// --8<-- [start:PermutationGroup]
// Define a cyclic group of order 3
Permutation p1 = {0, 1, 2};
Permutation p2 = {1, 2, 0};
Permutation p3 = {2, 0, 1};
auto C3 = PermutationGroup({p1, p2, p3});

XDIAG_PRINT(C3.size());
XDIAG_PRINT(C3.n_sites());
XDIAG_PRINT(C3.inverse(1)); // = 2
// --8<-- [end:PermutationGroup]
}


{
// --8<-- [start:Representation]
Representation r1 = {1, -1, 1, -1};
Representation r2 = {1, 1i, -1, -1i};

XDIAG_PRINT(r1 * r2);
// --8<-- [end:Representation]
}

  // clang-format on

} catch (Error e) {
  error_trace(e);
}
