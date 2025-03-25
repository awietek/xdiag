#include <vector>
#include <complex>
#include <math.h>
#include <xdiag/all.hpp>

using namespace xdiag;


std::complex<double> C_N_character(int N, int k, int p){
    return std::exp( ( (2 * k * p) * M_PI / (double)N ) * std::complex<double>(0,1.) );
};


int main(){
    try{

        // specify system parameters
        int N = 18;
        double J = 1.0;

        // fix number of up spins
        int Nup = N/2;

        // build Hamiltonian
        auto H = OpSum();
        for (int i=0; i<N; i++){
            H += "J" * Op("SdotS", {i, (i+1)%N});
        }
        H["J"] = J;

        // set up translation symmetry 
        // step 1: define shift by one site
        std::vector<int64_t> T_perm (N, 0);
        for (int i=0; i<N-1; i++){
            T_perm[i] = i + 1;
        };
        auto T = Permutation(T_perm);

        // step 2: define cyclic group
        std::vector<Permutation> perms (N);
        for (int i=0; i<N; i++){
            perms[i] = pow(T, i);
        };
        auto C_N = PermutationGroup(perms); 
        
        // step 3: define irreps from character tables, labelled by momentum (2pi/N ×) k
        std::vector<Representation> irreps (N);
        std::vector<std::complex <double> > characters (N, 0.0);
        for (int k = 0; k < N; k++){
            for (int p=0; p<N; p++){
                characters[p] = C_N_character(N, k, p);
            };
            irreps[k] =  Representation(C_N, characters);
        };
       
        // sort states in Sz_tot = 0 sector by translation irrep
        int lvl_per_block = 5;
        std::vector<double> TOS_x_vals(N*lvl_per_block); 
        std::vector<double> TOS_y_vals(N*lvl_per_block); 
        for (int k = 0; k < N; k++){
            auto sym_block = Spinhalf(N, Nup, irreps[k]);
            // set neigvals = ... for reasonable convergence of lowest eigenvalues
            int neigvals = 3;
            auto lanczos_result = eigvals_lanczos(H, sym_block, neigvals);
            for (int i=0; i<lvl_per_block; i++){
                TOS_x_vals[k*lvl_per_block + i] = (2 * M_PI * k) / N;
                TOS_y_vals[k*lvl_per_block + i] = lanczos_result.eigenvalues[i];
            }
        };

        // use the vectors TOS_x_vals and TOS_y_vals to create e.g. the TOS plot ...
        // ... for a simple plotting method check the Julia version of this example
        
        return 0;

    }catch(Error e){
        error_trace(e);
    }

}
