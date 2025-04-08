#include <xdiag/all.hpp>

using namespace xdiag;


std::complex<double> C_N_character(int N, int k, int p){
    return std::exp( ( (2 * k * p) * M_PI / (double)N ) * std::complex<double>(0,1.) );
};

std::vector<double> compute_level_statistics(int N, OpSum H){
    
    // fix magnetization (Sztot = 0 still has spin-flip symmetry!)
    int Nup = N/2 + 1;

    // fix lattice momentum (k = 0, N/2 still have parity symmetry!)
    int k = 1;

    std::vector<int64_t> T_perm (N, 0);
    for (int i=0; i<N-1; i++){
        T_perm[i] = i + 1;
    };
    auto T = Permutation(T_perm);

    std::vector<Permutation> perms (N);
    for (int i=0; i<N; i++){
        perms[i] = pow(T, i);
    };
    auto C_N = PermutationGroup(perms);

    std::vector< std::complex <double> > irrep_k_characters (N, 0.0);
    for (int p=0; p<N; p++){
        irrep_k_characters[p] = C_N_character(N, k, p);
    };
    auto irrep_k = Representation(C_N, irrep_k_characters);

    // block of Hamiltonian without remaining symmetries
    auto block = Spinhalf(N, Nup, irrep_k);

    // find its eigenspectrum
    arma::cx_mat Hmat = matrixC(H, block);     
    arma::vec eigenvalues = eig_sym(Hmat);

    // compute level statistics (taking only the inner most half of spectrum)
    int N_levels = eigenvalues.n_elem;
    int s_start = (int)(0.25 * N_levels);
    int s_stop = (int)(0.75 * N_levels);
    int s_num = s_stop - s_start;
    std::vector<double> s_arr (s_num);
    double s_arr_sum = 0.0;
    for (int i=0; i<s_num; i++){
        s_arr[i] = eigenvalues[s_start+i+1] - eigenvalues[s_start+i];
        s_arr_sum += s_arr[i];
    }
    s_arr_sum = s_arr_sum / s_num;

    // normalize
    for (int i=0; i<s_num; i++){
        s_arr[i] = s_arr[i] /  s_arr_sum;
    }

    return s_arr;
}


int main() try {
    
    int N = 18; // length of spin chain

    // definition of integrable model
    auto H_i = OpSum();
    for (int i=0; i<N; i++){
        H_i += "J" * Op("SzSz", {i, (i+1)%N});
        H_i += "Delta" * Op("Exchange", {i, (i+1)%N});
    }

    // definition of non-integrable model
    auto H_ni = OpSum();
    for (int i=0; i<N; i++){
        H_ni += "J" * Op("SzSz", {i, (i+1)%N});
        H_ni += "Delta" * Op("Exchange", {i, (i+1)%N});
        H_ni += "J2" * Op("SzSz", {i, (i+2)%N});
    }

    // assign coupling values
    double J, Delta, J2;
    J = 1.0;
    Delta = J2 = 0.5;

    H_i["J"] = H_ni["J"] = J;
    H_i["Delta"] = H_ni["Delta"] = Delta;
    H_ni["J2"] = J2;

    // compute level statistics (remember to eliminate trivial symmetries!)
    std::vector<double> H_i_statistics = compute_level_statistics(N, H_i); 
    std::vector<double> H_ni_statistics = compute_level_statistics(N, H_ni); 

    // do something with H_i_statistics and H_ni_statistics here (Julia version has a plotting routine)
    return 0;

} catch(Error e) {
    error_trace(e);
}


