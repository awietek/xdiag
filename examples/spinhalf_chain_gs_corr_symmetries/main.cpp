#include <xdiag/all.hpp>
#include <vector>
#include <complex>
#include <math.h>

using namespace xdiag;

// treat Hilbert space as a whole (memory consuming)
void gs_correlator_simple(int N, OpSum& H, OpSum& corr_op){
    auto full_block = Spinhalf(N);
    auto lanczos_res = eigs_lanczos(H, full_block);
    auto e0 = lanczos_res.eigenvalues[0];
    auto psi0 = lanczos_res.eigenvectors;
    auto corr = inner(corr_op, psi0);
   
    Log("Ground state energy: {:.12f}", e0);
    Log("Ground state correlator: {:.12f}", corr);
    Log("State vector length: {:d}", size(psi0));
};

// use Sz_tot conservation (less memory consuming)
void gs_correlator_Sz_sym(int N, OpSum& H, OpSum& corr_op){
    double e0 = 0.0;
    double gs;
    State psi0;
    // check blocks with fixed number of up-spins
    for (int nup = 0; nup < N+1; nup++){
        auto sym_block = Spinhalf(N, nup);
        auto lanczos_res = eigs_lanczos(H, sym_block);
        gs = lanczos_res.eigenvalues[0];
        if (gs < e0){
            e0 = gs;
            psi0 = lanczos_res.eigenvectors;
        };
    };
    auto corr = inner(corr_op, psi0);

    Log("Ground state energy: {:.12f}", e0);
    Log("Ground state correlator: {:.12f}", corr);
    Log("State vector length: {:d}", size(psi0));
};

std::complex<double> C_N_character(int N, int k, int p){
    return std::exp( ( (2 * k * p) * M_PI / (double)N ) * std::complex<double>(0,1.) );

};

// use translation symmetry (less memory consuming)
void gs_correlator_translation_sym(int N, OpSum& H, OpSum& corr_op){

    // define shift by one site
    std::vector<int64_t> T_perm (N, 0);
    for (int i=0; i<N-1; i++){
        T_perm[i] = i + 1;
    };
    auto T = Permutation(T_perm);

    // define cyclic group
    std::vector<Permutation> perms (N);
    for (int i=0; i<N; i++){
        perms[i] = pow(T, i);
    };
    auto C_N = PermutationGroup(perms);
    
    
    // define irreps from character tables, labelled by momentum (2pi/N ×) k
    std::vector<Representation> irreps (N);
    std::vector<std::complex <double> > characters (N, 0.0);
    for (int k = 0; k < N; k++){
        for (int p=0; p<N; p++){
            characters[p] = C_N_character(N, k, p);
        };
        irreps[k] =  Representation(C_N, characters);
    };
   
    // check all translation symmetry blocks
    double e0 = 0.0;
    double gs;
    State psi0;
    for (int k = 0; k < N; k++){
        auto sym_block = Spinhalf(N, irreps[k]);
        auto lanczos_res = eigs_lanczos(H, sym_block);
        gs = lanczos_res.eigenvalues[0];
        if (gs < e0){
            e0 = gs;
            psi0 = lanczos_res.eigenvectors;
        };
    };

    // now the operator must be symmetrized!
    auto corr = inner(symmetrize(corr_op, C_N), psi0);
    
    Log("Ground state energy: {:.12f}", e0);
    Log("Ground state correlator: {:.12f}", corr);
    Log("State vector length: {:d}", size(psi0));
};




int main() try {

    // periodic N-site Heisenberg antiferromagnet
    int N = 14;
    auto H = OpSum();
    for (int i = 0; i < N; i++) {
        H += "J" * Op("SdotS", {i, (i+1) % N});
    };
    H["J"] = 1.0;

    // define spin correlator at "half-chain" distance
    auto corr_op = OpSum();
    corr_op = Op("SdotS", {0, N/2 - 1});

    // run different implementations
    gs_correlator_simple(N, H, corr_op);
    
    gs_correlator_Sz_sym(N, H, corr_op);

    gs_correlator_translation_sym(N, H, corr_op);
    
} catch (Error e) {
        error_trace(e);
};
