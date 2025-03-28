#include <xdiag/all.hpp>

using namespace xdiag;

int main()
{
    say_hello();
    unsigned Lx = 4;
    unsigned Ly = 4;

    // construct Hamiltonian
    auto ops = OpSum();
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
    ops += "U" * Op("HubbardU");

    ops["t"] = -1.0;
    ops["U"] = 0.0;

    unsigned nup = Lx*Ly/2;
    unsigned ndn = Lx*Ly/2;
    auto block = Electron(Lx*Ly, nup, ndn);

    auto [e0,psi] = eig0(ops,block);

    ops["U"] = -10.0;

    auto obs = OpSum();
    obs += Op("HubbardU");

    arma::vec times = arma::linspace(0,100,2);
    arma::vec obs_vec = arma::vec(times.size());

    for (unsigned i=0; i<times.size(); i++)
    {
        time_evolve_inplace(ops,psi,times[1]-times[0]);
        obs_vec[i] = inner(obs,psi);
    }

    // save and plot

    return 0;
}