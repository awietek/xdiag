#include <xdiag/all.hpp>

using namespace xdiag;


int main() 
try {

    // define open ferromagnetic XXZ chain
    int N = 16;
    double J = 0.1;
    double Delta = 0.5;

    auto H = OpSum();
    for (int i=0; i<(N-1); i++){
        H += "J" * Op("SzSz", {i, i+1});
        H += "Delta" * Op("Exchange", {i, i+1});
    }

    H["J"] = J;
    H["Delta"] = Delta;

    // define initial state with domain wall
    auto block = Spinhalf(N);
    std::vector<std::string> psi0_vec (N);
    for (int i=0; i<N/2; i++) {psi0_vec[i] = "Up";}
    for (int i=N/2; i<N; i++) {psi0_vec[i] = "Dn";}
    auto psi = product_state(block, psi0_vec);

    // time evolve psi0 and measure Sz expectation value
    double dt = 0.5;
    int Nt = 30;
    std::vector<std::vector<double> > Sz_expectation(Nt, std::vector<double>(N));
    for (int t_step=0; t_step<Nt; t_step++){
        time_evolve_inplace(H, psi, dt);
        for (int i=0; i<N; i++){
            Sz_expectation[t_step][i] = innerC(Op("Sz", {i}), psi).real(); // result must be real
        }
    }

    // do something here with Sz_expectation (see Julia version for plotting routine)

    return 0;
} catch(Error e){
    error_trace(e);
}
