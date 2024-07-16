
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
 
  // clang-format on

} catch (Error e) {
  error_trace(e);
}
