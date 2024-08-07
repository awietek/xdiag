
#include <xdiag/all.hpp>

using namespace xdiag;

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
XDIAG_SHOW(C3.n_sites());
XDIAG_SHOW(C3.inverse(1)); // = 2
// --8<-- [end:PermutationGroup]
}


{
// --8<-- [start:Representation]
Representation r1 = {1, -1, 1, -1};
Representation r2 = {1, 1i, -1, -1i};

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
auto irrep = Representation({1, -1, 1, -1});
auto block_sym = Spinhalf(N, group, irrep);
XDIAG_SHOW(block_sym);

// with symmetries and Sz
auto block_sym_sz = Spinhalf(N, nup, group, irrep);
XDIAG_SHOW(block_sym_sz);

XDIAG_SHOW(block_sym_sz.n_sites());
XDIAG_SHOW(block_sym_sz.size());

// Iteration
for (auto pstate : block_sym_sz) {
  Log("{} {}", to_string(pstate), block_sym_sz.index(pstate));
}
XDIAG_SHOW(block_sym_sz.permutation_group());
XDIAG_SHOW(block_sym_sz.irrep());
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
auto irrep = Representation({1, -1, 1, -1});
auto block_sym = tJ(N, nup, ndn, group, irrep);
XDIAG_SHOW(block_sym);

XDIAG_SHOW(block_sym.n_sites());
XDIAG_SHOW(block_sym.size());

// Iteration
for (auto pstate : block_sym) {
  Log("{} {}", to_string(pstate), block_sym.index(pstate));
}
XDIAG_SHOW(block_sym.permutation_group());
XDIAG_SHOW(block_sym.irrep());
// --8<-- [end:tJ]
}


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
auto irrep = Representation({1, -1, 1, -1});
auto block_sym = Electron(N, group, irrep);
XDIAG_SHOW(block_sym);

// with symmetries and number conservation
auto block_sym_np = Electron(N, nup, ndn, group, irrep);
XDIAG_SHOW(block_sym_np);

XDIAG_SHOW(block_sym_np.n_sites());
XDIAG_SHOW(block_sym_np.size());

// Iteration
for (auto pstate : block_sym_np) {
  Log("{} {}", to_string(pstate), block_sym_np.index(pstate));
}
XDIAG_SHOW(block_sym_np.permutation_group());
XDIAG_SHOW(block_sym_np.irrep());
// --8<-- [end:Electron]

{
// --8<-- [start:matrix]
// Creates matrix H_{k=2} in Eq (18.23) of https://link.springer.com/content/pdf/10.1007/978-3-540-74686-7_18.pdf
int N = 4;
int nup = 3;
int ndn = 2;

// Define a Hubbard chain model
auto ops = OpSum();
for (int i=0; i< N; ++i){
  ops += Op("HOP", "T", {i, (i+1) % N});
}
ops["T"] = 1.0;
ops["U"] = 5.0;

// Create the a permutation group
auto p1 = Permutation({0, 1, 2, 3});
auto p2 = Permutation({1, 2, 3, 0});
auto p3 = Permutation({2, 3, 0, 1});
auto p4 = Permutation({3, 0, 1, 2});
auto group = PermutationGroup({p1, p2, p3, p4});
auto irrep = Representation({1, -1, 1, -1});
auto block = Electron(N, nup, ndn, group, irrep);

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
  ops += Op("HB", "J", {i, (i+1) % N});
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
  ops += Op("HB", "J", {i, (i+1) % N});
}
ops["J"] = 1.0;
auto [e0, gs] = eig0(ops, block);
// --8<-- [end:eig0]
}



{
// --8<-- [start:op]
auto op = Op("HOP", "T", {1, 2});
XDIAG_SHOW(op);
XDIAG_SHOW(op.type());
XDIAG_SHOW(op.coupling().as<std::string>());
XDIAG_SHOW(op.size());
XDIAG_SHOW(op[0]);
XDIAG_SHOW(op[1]);
XDIAG_SHOW(op.isexplicit());

 op = Op("HOP", 1.23, {1, 2});
XDIAG_SHOW(op);
XDIAG_SHOW(op.isreal());
XDIAG_SHOW(op.ismatrix());
XDIAG_SHOW(op.isexplicit());

arma::cx_mat m(arma::mat("0 0; 0 0"), arma::mat("0 -1; 1 0"));
op = Op("SY", m, 1);
XDIAG_SHOW(op);
XDIAG_SHOW(op.isreal());
XDIAG_SHOW(op.ismatrix());
XDIAG_SHOW(op.isexplicit());
// --8<-- [end:op]


// --8<-- [start:opsum]
// Define the 1D transverse-field Ising chain
int N = 12;
double J = 1.0;
double h = 0.5;
auto Sx = arma::mat("0 1; 1 0");

auto ops = OpSum();
for (int i = 0; i<N; ++i) {
  ops += Op("ISING", "J", {i, (i+1)%N});
  ops += Op("SX", arma::mat(h*Sx), i);
}
ops["J"] = 1.0;
XDIAG_SHOW(ops);
XDIAG_SHOW(ops.defined("J"));
XDIAG_SHOW(ops.isreal());
XDIAG_SHOW(ops.isexplicit());
// --8<-- [end:opsum]
}



{
// --8<-- [start:coupling]
auto cpl = Coupling("J");
XDIAG_SHOW(cpl.type());
XDIAG_SHOW(cpl.isexplicit());

cpl = Coupling(1.23);
XDIAG_SHOW(cpl.ismatrix());
XDIAG_SHOW(cpl.as<double>());
XDIAG_SHOW(cpl.as<complex>());

cpl = Coupling(arma::mat("1 2; -2 1"));
XDIAG_SHOW(cpl.ismatrix());
XDIAG_SHOW(cpl.isreal());
XDIAG_SHOW(cpl.as<arma::mat>());
XDIAG_SHOW(cpl.as<arma::cx_mat>());
// --8<-- [end:coupling]
} 



{
// --8<-- [start:state]
auto block = Spinhalf(2);
auto psi1 = State(block, arma::vec("1.0 2.0 3.0 4.0"));
XDIAG_SHOW(psi1);
XDIAG_SHOW(psi1.vector());
psi1.make_complex();
XDIAG_SHOW(psi1.vectorC());

auto psi2 = State(block, false, 3);
XDIAG_SHOW(psi2);
XDIAG_SHOW(psi2.matrixC());

auto psi3 = State(block, arma::cx_vec(arma::vec("1.0 2.0 3.0 4.0"),
				      arma::vec("4.0 3.0 2.0 1.0")));
XDIAG_SHOW(psi3.vectorC());
XDIAG_SHOW(psi3.real().vector());
XDIAG_SHOW(psi3.imag().vector());
// --8<-- [end:state]


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
XDIAG_SHOW(pstate.n_sites());
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
auto state = product(block, {"Up", "Dn"});
XDIAG_SHOW(state.vector());

zero(state);
XDIAG_SHOW(state.vector());

state = rand(block, false, 1234, true);
XDIAG_SHOW(state.vectorC());

state = zeros(block, true, 2);
XDIAG_SHOW(state.vector());
// --8<-- [end:create_state]
}


{
// --8<-- [start:algebra]
int N = 8;
auto block = Spinhalf(N,  N / 2);
auto ops = OpSum();
for (int i=0; i<N; ++i) {
  ops += Op("HB", 1.0, {i, (i+1)%N});
}
auto [e0, psi] = eig0(ops, block);

XDIAG_SHOW(norm(psi));
XDIAG_SHOW(norm1(psi));
XDIAG_SHOW(norminf(psi));

XDIAG_SHOW(dot(psi, psi));
XDIAG_SHOW(e0);
XDIAG_SHOW(inner(ops, psi));

auto phi = rand(block);
XDIAG_SHOW(phi.vector());
XDIAG_SHOW(psi.vector());
XDIAG_SHOW((psi + 2.0*phi).vector());
XDIAG_SHOW((psi*complex(0,3.0) + phi/2.0).vectorC());
// --8<-- [end:algebra]
}
 
 // clang-format on

} catch (Error e) {
  error_trace(e);
}
