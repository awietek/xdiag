#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {
  // clang-format off
 
{
// --8<-- [start:usage_guide_hs1]
int N = 8;
auto hs = Spinhalf(N);
// --8<-- [end:usage_guide_hs1]

// --8<-- [start:usage_guide_hs2]
for (auto spins : hs) {
  Log("{}", to_string(spins));
  Log("{}", index(hs, spins));
}
Log("dim: {}", size(hs));
// --8<-- [end:usage_guide_hs2]

// --8<-- [start:usage_guide_hs3]
int nup = 2;
auto b1 = Spinhalf(N, nup);

int ndn = 1;
auto b2 = tJ(N, nup, ndn);
auto b3 = Electron(N, nup, ndn);
// --8<-- [end:usage_guide_hs3]

// --8<-- [start:usage_guide_op1]
auto ops = OpSum();
for (int i=0; i<N; ++i) {
    int s1 = i;
    int s2 = (i+1) % N
    ops += "J" * Op("SdotS", {s1, s2});
}
ops["J"] = 1.0;
// --8<-- [end:usage_guide_op1]

// --8<-- [start:usage_guide_mat1]
arma::mat H = matrix(ops, block);
H.print();
// --8<-- [end:usage_guide_mat1]

// --8<-- [start:usage_guide_mat2]
arma::vec evals;
arma::mat evecs;
arma::eig_sym(evals, evecs, H);
// --8<-- [end:usage_guide_mat2]

// --8<-- [start:usage_guide_stat1]
bool real = true;
auto psi1 = State(b, real);
auto psi2 = zero_state(b, real);
// --8<-- [end:usage_guide_stat1]

// --8<-- [start:usage_guide_stat2]
int d = size(block)
arma::vec v(d, arma::fill::randu);
auto psi = State(b, v);
// --8<-- [end:usage_guide_stat2]

// --8<-- [start:usage_guide_stat3]
auto psi1 = product_state(block, {"Up", "Dn"});
auto psi2 = random_state(block);
// --8<-- [end:usage_guide_stat3]

// --8<-- [start:usage_guide_stat4]
double nrm = norm(psi);
double d = dot(psi1, psi2);
complex dc = dotC(psi1, psi2);
// --8<-- [end:usage_guide_stat4]

// --8<-- [start:usage_guide_stat5]
arma::vec v = vector(psi);
arma::cx_vec vc = vectorC(psi);
// --8<-- [end:usage_guide_stat5]

// --8<-- [start:usage_guide_stat6]
auto phi = apply(H, psi);
// --8<-- [end:usage_guide_stat6]

// --8<-- [start:usage_guide_iter1]
double e0 = eigval0(H, block);
// --8<-- [end:usage_guide_iter1]

// --8<-- [start:usage_guide_iter2]
auto [e0 , psi0] = eig0(H, block);
// --8<-- [end:usage_guide_iter2]

// --8<-- [start:usage_guide_iter3]
double t = 1.0;
auto phi = time_evolve(H, psi0, t);
// --8<-- [end:usage_guide_iter3]

// --8<-- [start:usage_guide_measu1]
for (int i=0; i<N; ++i) {
  auto op = Op("SzSz", {0, i});
  double corr = inner(op, psi0);
}
// --8<-- [end:usage_guide_measu1]

// --8<-- [start:usage_guide_io1]
auto fl = FileToml("spinhalf_chain.toml");
auto ops = read_opsum(fl, "Interactions");
ops["J"] = 1.0
// --8<-- [end:usage_guide_io1]

// --8<-- [start:usage_guide_io2]
auto f1 = FileH5("output.h5", "w!");
f1["e0"] = e0;
f1["evals"] = evals;
f1["evecs"] = evecs;
// --8<-- [end:usage_guide_io2]

// --8<-- [start:usage_guide_sym1]
auto T = Permutation({1, 2, 3, 4, 5, 6, 7, 0});
// --8<-- [end:usage_guide_sym1]

// --8<-- [start:usage_guide_sym2]
auto group = PermutationGroup({
  pow(T, 0), pow(T, 1), pow(T, 2), pow(T, 3),
  pow(T, 4), pow(T, 5), pow(T, 6), pow(T, 7)});
// --8<-- [end:usage_guide_sym2]

// --8<-- [start:usage_guide_sym3]
auto chi = arma::vec({1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0});
auto k = Representation(group, chi);
// --8<-- [end:usage_guide_sym3]

// --8<-- [start:usage_guide_sym4]
auto blk = Spinhalf(N, nup, irrep);
// --8<-- [end:usage_guide_sym4]

// --8<-- [start:usage_guide_sym5]
for (auto spins: blk) {
  Log("{}", to_string(spins));
}
// --8<-- [end:usage_guide_sym5]

// --8<-- [start:usage_guide_sym6]
auto fl = FileToml("symmetries.toml");
auto group = read_permutation_group(fl, "Symmetries");
auto irrep = read_representation(fl, "k.zero", "Symmetries");
// --8<-- [end:usage_guide_sym6]

// --8<-- [start:usage_guide_sym7]
auto og = symmetrize(ops, group);
auto oi = symmetrize(ops, irrep);
// --8<-- [end:usage_guide_sym7]

// --8<-- [start:usage_guide_dist1]
#include <xdiag/all.hpp>
using namespace xdiag;
int main(int argc, char* argv[]) try {
    MPI_Init(argc, argv);
    // genuine XDiag code here
    MPI_Finalize();
} catch (Error e) {
    error_trace(e);
}
// --8<-- [end:usage_guide_dist1]

// --8<-- [start:usage_guide_dist2]
auto block = SpinhalfDistributed(N, nup);
OpSum ops;    
for (int i = 0; i < N; ++i) {
    ops += Op("SdotS", {i, (i + 1) % N});
}
double e0 = eigval0(ops, block);
// --8<-- [end:usage_guide_dist2]

}

// clang-format on

} catch (Error e) {
  error_trace(e);
}
