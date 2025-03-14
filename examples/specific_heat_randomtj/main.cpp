#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <xdiag/all.hpp>

const unsigned Nsites = 8;   // total number of sites
const unsigned Nup = 2;       // total number of spin-up fermions
const unsigned Ndn = 2;       // total number of spin-down fermions
const unsigned Nsamples = 200; // number of disorder averages
const double thopp = 1.0;     // Variances hopping elements
const double Jhopp = 1.0;     // Variances exchange elements

// Create random devices to generate random couplings
std::random_device rd;
std::vector<unsigned int> seeds = {rd(), rd(), rd(), rd(), rd()};
std::seed_seq seq(seeds.begin(), seeds.end());
std::mt19937 gen(seq);

using namespace xdiag;

void get_H(xdiag::OpSum &H);

int main()
try
{
    // define Hamiltonian
    auto block = tJ(Nsites, Nup, Ndn); // tj model sites

    // Temperature range to compute the specific heat
    arma::vec Temp = arma::linspace<arma::vec>(0.01, 0.5, 64);

    // array to store specific heat
    arma::mat C = arma::mat(Temp.n_elem, Nsamples, arma::fill::zeros);

    // Obtain different disorder realizations
    arma::vec eigs;
    arma::cx_mat vecs;

    // Generate bond interactions. Save all possible elements in a vector of strings to modifie latter
    auto H = OpSum();
    for (unsigned i = 0; i < Nsites; i++)
    {
        for (unsigned j = (i + 1); j < Nsites; j++)
        {
            std::string Sij = "S_" + std::to_string(i) + "_" + std::to_string(j);
            std::string tij = "t_" + std::to_string(i) + "_" + std::to_string(j);
            H += Sij * Op("SdotS", {i, j});
            H += tij * Op("Hop", {i, j});
        }
    }

    for (unsigned i = 0; i < Nsamples; i++)
    {

        get_H(H); // generate new random couplings for the tj Hamiltonian

        arma::cx_mat Hmat = matrixC(H, block); // convert to matrix
        arma::eig_sym(eigs, vecs, Hmat);       // perform exact diagonalization
        for (unsigned j = 0; j < Temp.n_elem; j++)
        {
            arma::vec exp_eval = arma::exp(-eigs / Temp(j));
            double Z = arma::sum(exp_eval); // Partition function;
            double energy = arma::sum(eigs % exp_eval) / Z;
            C(j, i) = (1.0 / pow(Temp(j), 2)) * (arma::sum(eigs % (eigs % exp_eval)) / Z - pow(energy, 2)); // Specific Heat
        }
    }

    std::string filename = "randomtj.Nsites." + std::to_string((unsigned)Nsites) + ".outfile.h5";

    auto save_fl = FileH5(filename, "w!");
    save_fl["Temperature"] = Temp;
    save_fl["Specific_Heat"] = C;

    return 0;
}
catch (Error e)
{
    error_trace(e);
}

void get_H(xdiag::OpSum &H)
{
    // Normal random generators
    std::normal_distribution<double> Jdist(0.0, pow(Jhopp, 2));
    std::normal_distribution<double> tdist(0.0, pow(thopp, 2)/ sqrt(2.0));

    // update random couplings
    for (unsigned i = 0; i < Nsites; i++)
    {
        for (unsigned j = (i + 1); j < Nsites; j++)
        {
            std::string Sij = "S_" + std::to_string(i) + "_" + std::to_string(j);
            std::string tij = "t_" + std::to_string(i) + "_" + std::to_string(j);
            H[Sij] = Jdist(gen) / sqrt(Nsites);
            H[tij] = std::complex<double>(tdist(gen), tdist(gen)) / sqrt(Nsites);
        }
    }
}