#include <iostream>
#include <string>
#include <xdiag/all.hpp>
#include <vector>
#include <algorithm>
#include <cstdio>
//#include <format>

using namespace xdiag;

int main() try 
{
    say_hello();
    unsigned Lx = 4;
    unsigned Ly = 4;
    unsigned Nmin = 0;
    unsigned Nmax = 2*Lx*Ly;
    std::vector<double> Us = {0.0,-2.0,-10.0};
    
    for (unsigned i=0; i<Us.size(); i++)
    {
        double U = Us[i];
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
                ops += "mu" * Op("Ntot",s); // chemical potential
            }
        }
        ops["t"] = 1;
        ops["mu"] = -U/2;
        ops["U"] = U;

        for (unsigned N=Nmin; N<=Nmax; N+=2)
        {
            unsigned nup = N/2;
            unsigned ndn = N/2;
            auto block = Electron(Lx*Ly, nup, ndn); // create Hilbert space

            auto res = eigs_lanczos(ops,block);

            //std::string filename = std::format("data/tos_ahm/U({})_N({})_Lx({})_Ly({}).h5", U, N, Lx, Ly);
            char buffer[50];
            std::sprintf(buffer,"data/tos_ahm/U(%f)_N(%d)_Lx(%d)_Ly(%d).h5", U, N, Lx, Ly);
            auto f = FileH5(buffer,"w!");
            f["spectrum"] = res.eigenvalues;
        }

    }
    return 0;
}

catch (Error e) {
    error_trace(e);
}
