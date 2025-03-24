//#include <iostream>
//#include <string>
#include <xdiag/all.hpp>
//#include <vector>
//#include <algorithm>
#include <cstdio>

using namespace xdiag;

int main() try
{
    say_hello();
    double U = -10.0;
    unsigned Lx = 2;
    unsigned Ly = 2;

    // Construct Hamiltonian
        
    auto ops = OpSum();
        
    ops += "U" * Op("HubbardU"); // apply Hubbard interaction to all sites

    for (unsigned x=0; x<Lx; x++)
    {
        for (unsigned y=0; y<Ly; y++)
        {
            unsigned s = x*Ly + y;
            unsigned sx = ((x+1)%Lx)*Ly + y;
            unsigned sy = x*Ly + (y+1)%Ly;
            
            ops += "t" * Op("Hop",{s,sx}); // hopping in x-direction
            ops += "t" * Op("Hop",{s,sy}); // hopping in y-direction
        }
    }
    ops["t"] = 1;
    ops["U"] = U;

    // Create Hilbert space
    unsigned nup = Lx*Ly/2;
    unsigned ndn = Lx*Ly/2;
    auto block = Electron(Lx*Ly, nup, ndn);

    // Get ground state
    auto [e0,psi0] = eig0(ops,block);

    // build correlation matrix 

    arma::mat corr = arma::mat(Lx,Ly);
    for (unsigned x=0; x<Lx; x++)
    {
        for (unsigned y=0; y<Ly; y++)
        {
            unsigned s = x*Ly + y;
            if (s == 1)
            {
                continue;
            }
            auto op = Op("NtotNtot",{1,s});
            corr(s) = inner(op,psi0);
        }
    }

    //char buffer[50];
    //std::sprintf(buffer,"data/ahm_correlations/U(%f)_Lx(%d)_Ly(%d).h5", U, Lx, Ly);
    //auto f = FileH5(buffer,"w!");
    //f["correlator"] = corr;
    //f["psi"] = psi0;
    return 0;
}

catch (Error e) {
    error_trace(e);
}