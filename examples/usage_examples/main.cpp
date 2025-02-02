#include <xdiag/all.hpp>

using namespace xdiag;
using namespace std::complex_literals;

int main() try {
  // clang-format off

{
// --8<-- [start:Permutation]
Permutation p1 = {0, 2, 1, 3};
Permutation p2 = {2, 0, 1, 3};

XDIAG_SHOW(inverse(p1));
XDIAG_SHOW(p1*p2);
// --8<-- [end:Permutation]
}

{
// --8<-- [start:PermutationGroup]
// Define a cyclic group of order 3
Permutation p1 = {0, 1, 2};
Permutation p2 = {1, 2, 0};
Permutation p3 = {2, 0, 1};
auto C3 = PermutationGroup({p1, p2, p3});

XDIAG_SHOW(C3.size());
XDIAG_SHOW(C3.nsites());
XDIAG_SHOW(C3.inverse(1)); // = 2
// --8<-- [end:PermutationGroup]
}


{
// --8<-- [start:Representation]
Permutation p = {1, 2, 3, 0};
auto C4 = PermutationGroup({pow(p, 0), pow(p, 1), pow(p, 2), pow(p, 3)});
Representation r1(C4, arma::vec{1.0, -1.0, 1.0, -1.0});
Representation r2(C4, arma::cx_vec{1.0, 1.0i, -1.0, -1.0i});
XDIAG_SHOW(r1 * r2);
// --8<-- [end:Representation]
}


{
// --8<-- [start:Spinhalf]
int N = 4;
int nup = 2;

// without Sz conservation
auto block = Spinhalf(N);
XDIAG_SHOW(block);

// with Sz conservation
auto block_sz = Spinhalf(N, nup);
XDIAG_SHOW(block_sz);
 
// with symmetries, without Sz
Permutation p1 = {0, 1, 2, 3};
Permutation p2 = {1, 2, 3, 0};
Permutation p3 = {2, 3, 0, 1};
Permutation p4 = {3, 0, 1, 2};
auto group = PermutationGroup({p1, p2, p3, p4});
auto irrep = Representation(group, arma::vec({1, -1, 1, -1}));
auto block_sym = Spinhalf(N, irrep);
XDIAG_SHOW(block_sym);

// with symmetries and Sz
auto block_sym_sz = Spinhalf(N, nup, irrep);
XDIAG_SHOW(block_sym_sz);

XDIAG_SHOW(block_sym_sz.nsites());
XDIAG_SHOW(block_sym_sz.size());

// Iteration
for (auto pstate : block_sym_sz) {
  Log("{} {}", to_string(pstate), index(block_sym_sz, pstate));
}
// --8<-- [end:Spinhalf]
}

{
// --8<-- [start:tJ]
int N = 4;
int nup = 2;
int ndn = 1;

// without permutation symmetries
auto block = tJ(N, nup, ndn);
XDIAG_SHOW(block);

// with permutation symmetries
auto p1 = Permutation({0, 1, 2, 3});
auto p2 = Permutation({1, 2, 3, 0});
auto p3 = Permutation({2, 3, 0, 1});
auto p4 = Permutation({3, 0, 1, 2});
auto group = PermutationGroup({p1, p2, p3, p4});
auto irrep = Representation(group, arma::vec{1, -1, 1, -1});
auto block_sym = tJ(N, nup, ndn, irrep);
XDIAG_SHOW(block_sym);
XDIAG_SHOW(block_sym.nsites());
XDIAG_SHOW(block_sym.size());

// Iteration
for (auto pstate : block_sym) {
  Log("{} {}", to_string(pstate), index(block_sym, pstate));
}
// --8<-- [end:tJ]
}

{
// --8<-- [start:Electron]
int N = 4;
int nup = 2;
int ndn = 1;

// without number conservation
auto block = Electron(N);
XDIAG_SHOW(block);

// with number conservation
auto block_np = Electron(N, nup, ndn);
XDIAG_SHOW(block_np);

// with symmetries, without number conservation
auto p1 = Permutation({0, 1, 2, 3});
auto p2 = Permutation({1, 2, 3, 0});
auto p3 = Permutation({2, 3, 0, 1});
auto p4 = Permutation({3, 0, 1, 2});
auto group = PermutationGroup({p1, p2, p3, p4});
auto irrep = Representation(group, arma::vec{1, -1, 1, -1});
auto block_sym = Electron(N, irrep);
XDIAG_SHOW(block_sym);

// with symmetries and number conservation
auto block_sym_np = Electron(N, nup, ndn, irrep);
XDIAG_SHOW(block_sym_np);
XDIAG_SHOW(block_sym_np.nsites());
XDIAG_SHOW(block_sym_np.size());

// Iteration
for (auto pstate : block_sym_np) {
  Log("{} {}", to_string(pstate), index(block_sym_np, pstate));
}
// --8<-- [end:Electron]
}


{
// --8<-- [start:matrix]
// Creates matrix H_{k=2} in Eq (18.23) of https://link.springer.com/content/pdf/10.1007/978-3-540-74686-7_18.pdf
int N = 4;
int nup = 3;
int ndn = 2;

// Define a Hubbard chain model
auto ops = OpSum();
for (int i=0; i< N; ++i){
  ops += "T" * Op("Hop", {i, (i+1) % N});
}
ops+= "U" * Op("HubbardU");
ops["T"] = 1.0;
ops["U"] = 5.0;

// Create the a permutation group
auto p1 = Permutation({0, 1, 2, 3});
auto p2 = Permutation({1, 2, 3, 0});
auto p3 = Permutation({2, 3, 0, 1});
auto p4 = Permutation({3, 0, 1, 2});
auto group = PermutationGroup({p1, p2, p3, p4});
auto irrep = Representation(group, arma::vec{1.0, -1.0, 1.0, -1.0});
auto block = Electron(N, nup, ndn, irrep);
auto H = matrix(ops, block);
H.print();
// --8<-- [end:matrix]
}

{
// --8<-- [start:eigval0]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += "J" * Op("SdotS", {i, (i+1) % N});
}
ops["J"] = 1.0;
double e0 = eigval0(ops, block);
// --8<-- [end:eigval0]
}

{
// --8<-- [start:eig0]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += "J" * Op("SdotS", {i, (i+1) % N});
}
ops["J"] = 1.0;
auto [e0, gs] = eig0(ops, block);
// --8<-- [end:eig0]
}

{
// --8<-- [start:eigvals_lanczos]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += "J" * Op("SdotS", {i, (i+1) % N});
}
ops["J"] = 1.0;

// With random intial state
auto res = eigvals_lanczos(ops, block);
XDIAG_SHOW(res.alphas);
XDIAG_SHOW(res.betas);
XDIAG_SHOW(res.eigenvalues);

// With specific initial state
auto psi0 = product_state(block, {"Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"});
auto res2 = eigvals_lanczos(ops, psi0);
XDIAG_SHOW(res.alphas);
XDIAG_SHOW(res.betas);
XDIAG_SHOW(res.eigenvalues);
// --8<-- [end:eigvals_lanczos]
}
 
{
// --8<-- [start:eigs_lanczos]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += "J" * Op("SdotS", {i, (i+1) % N});
}
ops["J"] = 1.0;

// With random intial state
auto res = eigs_lanczos(ops, block);
XDIAG_SHOW(res.alphas);
XDIAG_SHOW(res.betas);
XDIAG_SHOW(res.eigenvalues);
XDIAG_SHOW(res.eigenvectors);
// --8<-- [end:eigs_lanczos]
}



{
// --8<-- [start:Op]
auto op = "T" * Op("Hop", {0, 1});
XDIAG_SHOW(op);

op = 1.23 * Op("Hop", {0, 1});
XDIAG_SHOW(op);

arma::cx_mat m(arma::mat("0 0; 0 0"), arma::mat("0 -1; 1 0"));
 op = Op("Matrix", 0, m);
XDIAG_SHOW(op);
XDIAG_SHOW(isreal(op));
// --8<-- [end:Op]
}

{
// --8<-- [start:OpSum]
// Define the 1D transverse-field Ising chain
int N = 12;
double J = 1.0;
double h = 0.5;
auto Sx = arma::mat("0 1; 1 0");

// Option 1: coupling constants as numbers
auto ops1 = OpSum();
for (int i = 0; i<N; ++i) {
  ops1 += J * Op("SzSz", {i, (i+1)%N});
  ops1 += h * Op("Matrix", i, Sx);
}

// Option 2: coupling constants as strings
auto ops2 = OpSum();
for (int i = 0; i<N; ++i) {
  ops2 += "J" * Op("SzSz", {i, (i+1)%N});
  ops2 += "h" * Op("Matrix", i, Sx);
}
ops2["J"] = J;
ops2["h"] = h;

XDIAG_SHOW(isapprox(ops1, ops2));
XDIAG_SHOW(isapprox(ops1 + ops2, 2.0 * ops1));
// --8<-- [end:OpSum]
}

{
// --8<-- [start:hc]
auto cdagup = Op("Cdagup", 0);
auto sdots = Op("SdotS", {0, 1});
auto hop = (1.0 + 1.0i) * Op("Hop", {0, 1});
XDIAG_SHOW(cdagup == hc(cdagup));
XDIAG_SHOW(sdots == hc(sdots));
XDIAG_SHOW(hop == hc(hop));
// --8<-- [end:hc]
}

{
// --8<-- [start:matrixtype1]
auto sx = arma::mat({{0, 1},{1, 0}});
auto op = Op("Matrix", 0, sx);
// --8<-- [end:matrixtype1]
}

{
// --8<-- [start:matrixtype2]
auto sx = arma::mat({{0, 1},{1, 0}});
auto sz = arma::mat({{0.5, 1},{0, -0.5}});

arma::mat sxsz = arma::kron(sx, sz);
arma::mat sxszsxsz = arma::kron(sxsz, sxsz);
 
auto op_sxsz = Op("Matrix", {0, 1}, sxsz);
auto op_sxszsxsz = Op("Matrix", {0, 1, 2, 3}, sxsz);
// --8<-- [end:matrixtype2]
}

{
// --8<-- [start:state]
auto block = Spinhalf(2);
auto psi1 = State(block, arma::vec("1.0 2.0 3.0 4.0"));
XDIAG_SHOW(psi1);
XDIAG_SHOW(vector(psi1));
make_complex(psi1);
XDIAG_SHOW(vectorC(psi1));

auto psi2 = State(block, false, 3);
XDIAG_SHOW(psi2);
XDIAG_SHOW(matrixC(psi2));

auto psi3 = State(block, arma::cx_vec(arma::vec("1.0 2.0 3.0 4.0"),
				      arma::vec("4.0 3.0 2.0 1.0")));
XDIAG_SHOW(vectorC(psi3));
XDIAG_SHOW(vector(real(psi3)));
XDIAG_SHOW(vector(imag(psi3)));
// --8<-- [end:state]
}

{
// --8<-- [start:product_state]
auto pstate = ProductState({"Up", "Dn", "Emp", "UpDn"});
for (auto s : pstate) {
  Log("{}", s);
}
XDIAG_SHOW(to_string(pstate));

pstate = ProductState();
pstate.push_back("Dn");
pstate.push_back("Up");
pstate.push_back("Dn");
XDIAG_SHOW(pstate.nsites());
for (auto s : pstate) {
  Log("{}", s);
}
XDIAG_SHOW(to_string(pstate));
// --8<-- [end:product_state]
}


{
// --8<-- [start:random_state]
auto block = Spinhalf(2);
auto state = State(block, false);  // complex State
auto rstate1 = RandomState(1234);
fill(state, rstate1);
XDIAG_SHOW(state.vectorC());

auto rstate2 = RandomState(4321);
fill(state, rstate2);
XDIAG_SHOW(state.vectorC());

fill(state, rstate1);
XDIAG_SHOW(state.vectorC());
// --8<-- [end:random_state]
}

{
// --8<-- [start:fill]
auto block = Spinhalf(2);
auto state = State(block);  
auto pstate = ProductState({"Up", "Dn"});
fill(state, pstate);
XDIAG_SHOW(state.vector());

auto rstate = RandomState(1234);
fill(state, rstate);
XDIAG_SHOW(state.vector());
// --8<-- [end:fill]
}


{
// --8<-- [start:create_state]
auto block = Spinhalf(2);
auto state = product_state(block, {"Up", "Dn"});
XDIAG_SHOW(state.vector());

zero(state);
XDIAG_SHOW(state.vector());

state = random_state(block, false, 1234, true);
XDIAG_SHOW(state.vectorC());

state = zero_state(block, true, 2);
XDIAG_SHOW(state.vector());
// --8<-- [end:create_state]
}

{
// --8<-- [start:time_evolve]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("SdotS", {i, (i+1) % N});
}

auto psi0 = product_state(block, {"Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"});
double time = 1.0;
auto psi = time_evolve(ops, psi0, time);
time_evolve_inplace(ops, psi0, time);
XDIAG_SHOW(isapprox(psi0, psi));
// --8<-- [end:time_evolve]
}

{
// --8<-- [start:imaginary_time_evolve]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("SdotS", {i, (i+1) % N});
}

// Compute ground state energy
double e0 = eigval0(ops, block);
 
auto psi0 = product_state(block, {"Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"});
double time = 1.0;
double precision = 1e-12;
auto psi = imaginary_time_evolve(ops, psi0, time, precision, e0);
imaginary_time_evolve_inplace(ops, psi0, time, precision, e0);
XDIAG_SHOW(isapprox(psi0, psi));
// --8<-- [end:imaginary_time_evolve]
}
 
 
{
// --8<-- [start:evolve_lanczos]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("SdotS", {i, (i+1) % N});
}

// Compute ground state energy
double e0 = eigval0(ops, block);
 
auto psi0 = product_state(block, {"Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"});
double time = 1.0;
double precision = 1e-12;
auto res = evolve_lanczos(ops, psi0, time, precision, e0, true, 500);
XDIAG_SHOW(res.alphas);
XDIAG_SHOW(res.betas);
// --8<-- [end:evolve_lanczos]
}
 
{
// --8<-- [start:time_evolve_expokit]
int N = 8;
int nup = N / 2;
auto block = Spinhalf(N, nup);
    
// Define the nearest-neighbor Heisenberg model
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("SdotS", {i, (i+1) % N});
}

auto psi0 = product_state(block, {"Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"});
double time = 1.0;
auto res1 = time_evolve_expokit(ops, psi0, time, 20);
auto res2 = time_evolve_expokit_inplace(ops, psi0, time, 20);
XDIAG_SHOW(isapprox(psi0, res1.state));
XDIAG_SHOW(res1.error);
XDIAG_SHOW(res1.hump);
// --8<-- [end:time_evolve_expokit]
}
 

 
{
// --8<-- [start:algebra]
int N = 8;
auto block = Spinhalf(N,  N / 2);
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("SdotS", {i, (i+1)%N});
}
auto [e0, psi] = eig0(ops, block);

XDIAG_SHOW(norm(psi));
XDIAG_SHOW(norm1(psi));
XDIAG_SHOW(norminf(psi));

XDIAG_SHOW(dot(psi, psi));
XDIAG_SHOW(e0);
XDIAG_SHOW(inner(ops, psi));

auto phi = random_state(block);
XDIAG_SHOW(phi.vector());
XDIAG_SHOW(psi.vector());
XDIAG_SHOW((psi + 2.0*phi).vector());
XDIAG_SHOW((psi*complex(0,3.0) + phi/2.0).vectorC());
// --8<-- [end:algebra]
}

{
// --8<-- [start:apply] 
int N = 8;
auto block = Spinhalf(N,  N / 2);
auto ops = OpSum();
for (int i=0; i<N; ++i){
  ops += Op("SdotS", {i, (i+1)%N});
}
auto [e0, psi] = eig0(ops, block);
auto phi = apply(Op("S+", 2), psi);
XDIAG_SHOW(inner(ops, psi));
XDIAG_SHOW(inner(ops, phi));
// --8<-- [end:apply]
}

{
// --8<-- [start:symmetrize]
int N = 4;
int nup = 2;
auto block = Spinhalf(N, nup);
auto p1 = Permutation({0, 1, 2, 3});
auto p2 = Permutation({1, 2, 3, 0});
auto p3 = Permutation({2, 3, 0, 1});
auto p4 = Permutation({3, 0, 1, 2});
auto group = PermutationGroup({p1, p2, p3, p4});
auto rep = Representation(group);
auto block_sym = Spinhalf(N, rep);

auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("SdotS", {i, (i+1)%N});
}
auto [e0, psi] = eig0(ops, block);
auto [e0s, psi_sym] = eig0(ops, block_sym);

auto corr = Op("SdotS", {0, 1});
auto nn_corr = inner(corr, psi);
auto corr_sym = symmetrize(corr, group);
auto nn_corr_sym = innerC(corr_sym, psi_sym);
XDIAG_SHOW(nn_corr);
XDIAG_SHOW(nn_corr_sym);
// --8<-- [end:symmetrize]
}

 
{
// --8<-- [start:FileToml]
auto fl = FileToml(XDIAG_DIRECTORY "/misc/data/toml/input.toml");
XDIAG_SHOW(defined(fl, "N"));

int N = fl["N"].as<int>();
int nup = fl["nup"].as<int>();
double J1 = fl["J1"].as<double>();
double J2 = fl["J2"].as<double>();
 
auto block = Spinhalf(N, nup);
auto H = OpSum();
for (int i=0; i<N; ++i){
  H += J1 * Op("SdotS", {i, (i+1)%N});
  H += J2 * Op("SdotS", {i, (i+2)%N});
}
double e0 = eigval0(H, block);
XDIAG_SHOW(e0);
// --8<-- [end:FileToml]
}

{
// --8<-- [start:read_opsum]
std::string file = XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
auto fl = FileToml(file);
auto ops = read_opsum(fl, "Interactions");
XDIAG_SHOW(ops);
// --8<-- [end:read_opsum]
}

{
// --8<-- [start:read_permutation_group]
std::string file = XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
auto fl = FileToml(file);
auto group = read_permutation_group(fl, "Symmetries");
XDIAG_SHOW(group);
// --8<-- [end:read_permutation_group]
}

 {
// --8<-- [start:read_representation]
std::string file = XDIAG_DIRECTORY "/misc/data/irreps.toml";
auto fl = FileToml(file);

auto k_0 = read_representation(fl, "k_0");
XDIAG_SHOW(k_0);
XDIAG_SHOW(isreal(k_0));

auto k_pi2 = read_representation(fl, "k_pi2");
XDIAG_SHOW(k_pi2);
XDIAG_SHOW(isreal(k_pi2));

auto k_pi = read_representation(fl, "k_pi");
XDIAG_SHOW(k_pi);
XDIAG_SHOW(isreal(k_pi));

auto k_pi2_half = read_representation(fl, "k_pi2_half");
XDIAG_SHOW(k_pi2_half);
XDIAG_SHOW(isreal(k_pi2_half));
// --8<-- [end:read_representation]
}
 
{
// --8<-- [start:FileH5]
std::string filename = XDIAG_DIRECTORY "/misc/data/hdf5/write.h5";
auto fl = FileH5(filename, "w!");

// Write output to the hdf5 file
fl["val"] = 12;
fl["test/to"] = 22;
fl["test/to2/group"] = 32;
fl["test/to3/group2/asdf"] = 42;

auto mat = arma::cx_mat(3, 5, arma::fill::randn);
fl["a/b/c/mat"] = mat;
// --8<-- [end:FileH5]
}
 
 
// }
  // clang-format on
} catch (Error e) {
  error_trace(e);
}
