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
XDIAG_SHOW(size(hspace));
XDIAG_SHOW(size(block));
// --8<-- [end:first_steps_4]

// --8<-- [start:first_steps_5]
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += "J" * Op("SdotS", {i, (i+1) % N});
}
ops["J"] = 1.0;
// --8<-- [end:first_steps_5]

// --8<-- [start:first_steps_6]
auto [e0, psi0] = eig0(ops, block);
Log("e0: {:.12f}", e0);
// --8<-- [end:first_steps_6]

// --8<-- [start:first_steps_7]
arma::mat H = matrix(ops, block); 
arma::vec evals;
arma::mat evecs;
arma::eig_sym(evals, evecs, H);
Log("e0: {:.12f}, e1: {:.12f}", evals[0], evals[1]);
// --8<-- [end:first_steps_7]

// --8<-- [start:first_steps_8]
for (int i=1; i<N; ++i) {
  auto op = Op("SzSz", {0, i});
  double corr = inner(op, psi0);
  Log("<Sz_0 Sz_{}> = {:.12f}", i, corr);
}
// --8<-- [end:first_steps_8]

 
}

// clang-format on

} catch (Error e) {
  error_trace(e);
}
