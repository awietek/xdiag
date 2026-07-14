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
  Log("{}", to_string(spins));      // local quantum numbers as integers
  Log("{}", to_string(spins, hs));  // human-readable configuration
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

int d = 4;         // local dimension (occupations 0, 1, 2, 3)
int nbosons = 3;
auto b4 = Boson(N, d, nbosons);

int nfermions = 3;
auto b5 = Fermion(N, nfermions);
// --8<-- [end:usage_guide_hs3]

// --8<-- [start:usage_guide_op1]
auto ops = OpSum();
for (int i=0; i<N; ++i) {
    int s1 = i;
    int s2 = (i+1) % N;
    ops += "J" * Op("SdotS", {s1, s2});
}
ops["J"] = 1.0;
// --8<-- [end:usage_guide_op1]

// --8<-- [start:usage_guide_op2]
// "J" is a named coupling. Its value is set with the [] operator, ...
ops["J"] = 1.0;
// ... and plain() substitutes every named coupling by its numerical value.
auto ops_plain = ops.plain();
// --8<-- [end:usage_guide_op2]

// --8<-- [start:usage_guide_op3]
// OpSums can be multiplied, forming the (non-commutative) operator algebra product
auto sz0 = 1.0 * Op("Sz", 0);
auto sz1 = 1.0 * Op("Sz", 1);
auto szsz = sz0 * sz1;
// --8<-- [end:usage_guide_op3]

// --8<-- [start:usage_guide_op4]
// hermitian conjugation via hc(), e.g. hc(S+) = S-
auto sp = 1.0 * Op("S+", 0);
auto sm = hc(sp);
auto sx = sp + hc(sp);   // a hermitian combination
// --8<-- [end:usage_guide_op4]

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
int64_t d = size(block);

arma::vec v(d, arma::fill::randu);       // real coefficients
auto psi = State(block, v);

arma::cx_vec vc(d, arma::fill::randu);   // complex coefficients
auto psic = State(block, vc);
// --8<-- [end:usage_guide_stat1]

// --8<-- [start:usage_guide_stat2]
bool real = true;
auto psi1 = State(block, real);        // a zero state (created implicitly)
auto psi2 = zero_state(block, real);   // a zero state (created explicitly)
// --8<-- [end:usage_guide_stat2]

// --8<-- [start:usage_guide_stat3]
auto phi = product_state(block, {1, 0});   // per-site local states: Up=1, Dn=0
// --8<-- [end:usage_guide_stat3]

// --8<-- [start:usage_guide_stat_rand]
int64_t seed = 1234;
auto chi = random_state(block, true, 1, seed);
// --8<-- [end:usage_guide_stat_rand]

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
double e0 = eigval0(ops, block);
// --8<-- [end:usage_guide_iter1]

// --8<-- [start:usage_guide_iter2]
auto [e0, psi0] = eig0(ops, block);
// --8<-- [end:usage_guide_iter2]

// --8<-- [start:usage_guide_iter_lanczos]
auto res = eigs_lanczos(ops, block);
arma::vec eigenvalues = res.eigenvalues;   // Ritz eigenvalue estimates
State eigenvectors = res.eigenvectors;
arma::vec alphas = res.alphas;             // diagonal of the tridiagonal matrix
arma::vec betas = res.betas;               // off-diagonal of the tridiagonal matrix
// --8<-- [end:usage_guide_iter_lanczos]

// --8<-- [start:usage_guide_iter_eigvals]
int64_t neigs = 3;
arma::vec eigenvalues = eigvals(ops, block, neigs);
// --8<-- [end:usage_guide_iter_eigvals]

// --8<-- [start:usage_guide_iter_eigs]
auto [eigenvalues, eigenvectors] = eigs(ops, block, neigs);
// --8<-- [end:usage_guide_iter_eigs]

// --8<-- [start:usage_guide_iter_lobpcg]
auto res = eigs_lobpcg(ops, block, neigs);
arma::vec eigenvalues = res.eigenvalues;    // the neigs lowest eigenvalues
State eigenvectors = res.eigenvectors;      // corresponding eigenvectors
arma::vec residuals = res.residual_norms;   // final residual norm of each
// --8<-- [end:usage_guide_iter_lobpcg]

// --8<-- [start:usage_guide_iter3]
double t = 1.0;
auto phi = time_evolve(ops, psi0, t);
// --8<-- [end:usage_guide_iter3]

// --8<-- [start:usage_guide_time_imag]
double tau = 1.0;
auto eta = imaginary_time_evolve(ops, psi0, tau);
// --8<-- [end:usage_guide_time_imag]

// --8<-- [start:usage_guide_measu1]
for (int i=0; i<N; ++i) {
  auto op = Op("SzSz", {0, i});
  double corr = inner(op, psi0);
}
// --8<-- [end:usage_guide_measu1]

// --8<-- [start:usage_guide_measu2]
// Build an arbitrary composite operator via the algebra product, e.g. S^x_0 S^y_1
auto op = Op("Sx", 0) * Op("Sy", 1);
complex corr = innerC(op, psi0);
// --8<-- [end:usage_guide_measu2]

// --8<-- [start:usage_guide_expect]
// <psi0| Sz_i |psi0> on every site i, returned as a vector
arma::vec sz = expect(psi0, "Sz");
// --8<-- [end:usage_guide_expect]

// --8<-- [start:usage_guide_corr]
// C(i, j) = <psi0| Sz_i Sz_j |psi0> for all pairs of sites
arma::mat szsz = correlation_matrix(psi0, "Sz", "Sz");
// --8<-- [end:usage_guide_corr]

// --8<-- [start:usage_guide_io1]
auto fl = FileToml("spinhalf_chain.toml");
auto ops = read_opsum(fl, "Interactions");
ops["J"] = 1.0;
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
auto blk = Spinhalf(N, nup, k);
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

// --8<-- [start:usage_guide_spm1]
auto coo_mat = coo_matrix(ops, block);
auto csr_mat = csr_matrix(ops, block);
auto csc_mat = csc_matrix(ops, block);
// --8<-- [end:usage_guide_spm1]

// --8<-- [start:usage_guide_spm2]
auto csc_mat = csc_matrix(ops, block);
auto colptr = arma::conv_to<arma::uvec>::from(csc_mat.colptr);
auto row = arma::conv_to<arma::uvec>::from(csc_mat.row);
auto A = arma::sp_mat(colptr, row, csc_mat.data, csc_mat.nrows, csc_mat.ncols);
// --8<-- [end:usage_guide_spm2]

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
