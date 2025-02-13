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

{
// --8<-- [start:io_1]
auto fl = FileToml(XDIAG_DIRECTORY "/examples/user_guide/spinhalf_chain.toml");
auto ops = read_opsum(fl, "Interactions");
// --8<-- [end:io_1]
}

// --8<-- [start:io_2]
auto fl = FileH5(XDIAG_DIRECTORY "/examples/user_guide/output.h5", "w!");
fl["e0"] = e0;
fl["evals"] = evals;
fl["evecs"] = evecs; 
// --8<-- [end:io_2]


// --8<-- [start:symmetries_1]
auto T = Permutation({1, 2, 3, 4, 5, 6, 7, 0});
// --8<-- [end:symmetries_1]


// --8<-- [start:symmetries_2]
auto group = PermutationGroup({pow(T, 0), pow(T, 1), pow(T, 2), pow(T, 3),
                               pow(T, 4), pow(T, 5), pow(T, 6), pow(T, 7)});
// --8<-- [end:symmetries_2]

// --8<-- [start:symmetries_3]
auto irrep_k_0 = Representation(group, arma::vec{1.0, 1.0, 1.0, 1.0,
						 1.0, 1.0, 1.0, 1.0});
auto irrep_k_pi = Representation(group, arma::vec{1.0, -1.0, 1.0, -1.0,
						  1.0, -1.0, 1.0, -1.0});
// --8<-- [end:symmetries_3]
 
// --8<-- [start:symmetries_4]
auto block_k_0 = Spinhalf(N, nup, irrep_k_0);
auto block_k_pi = Spinhalf(N, nup, irrep_k_pi);
double e0_k_0 = eigval0(ops, block_k_0);
double e0_k_pi = eigval0(ops, block_k_pi);
Log("e0: k=0: {:.12f}, k=pi: {:.12f}", e0_k_0, e0_k_pi);
// --8<-- [end:symmetries_4]

}

// clang-format on

} catch (Error e) {
  error_trace(e);
}
