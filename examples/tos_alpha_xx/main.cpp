#include <iostream>
#include <string>
#include <xdiag/all.hpp>
using namespace xdiag;

const unsigned Nsites = 14; // total number of sites
const double alpha = 1.0;   // Exponent power-law decay

int main()
{
    auto ops = OpSum(); // create OpSum

    for (int i = 0; i < Nsites; i++)
    {
        for (int j = i + 1; j < Nsites; j++)
        {
            double J = -1.0 / (sqrt(Nsites) * pow(abs((double)(i - j)), alpha));
            ops += J * Op("Exchange", {i, j});
        }
    }

    // Diagonalize in each magnetization sector using FullED
    unsigned total_eigen = pow(2, Nsites);
    arma::mat energies = arma::mat(total_eigen, 2, arma::fill::zeros); // Collect the energies with the number of spins up.
    arma::vec eig;

    unsigned count = 0;
    for (unsigned nup = 0; nup <= Nsites; nup++)
    {
        auto block = Spinhalf(Nsites, nup);
        arma::mat H = matrix(ops, block);
        eig_sym(eig, H);
        for (unsigned j = 0; j < eig.n_elem; j++)
        {
            energies(count, 0) = nup;    // Magnetization sector
            energies(count, 1) = eig(j); // Energy eigenvalue
            count++;
        }
    }
    // Construct the filename
    std::string flstring = "energies_tos_XXmodel.Nsites." + std::to_string(Nsites) + ".alpha" + std::to_string((double)alpha) + ".outfile.h5";

    auto save_fl = FileH5(flstring, "w!");
    save_fl["energies"] = energies;
    return 0;
}