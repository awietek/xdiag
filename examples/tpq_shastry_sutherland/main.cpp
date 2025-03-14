#include <iostream>
#include <vector>
#include <string>
#include <xdiag/all.hpp>

const unsigned Nsites = 20; // Total number of sites in the Shastry-Sutherland model:
const unsigned Rtqp = 5;    // Number of random vectors to be used in the algorithm

auto fl = xdiag::FileToml("shastry_sutherland_L_5_W_4.toml"); // TOML file with the list of interactions
using namespace xdiag;

void Tqp_mag_sector(arma::mat &Obser, arma::vec &DBetas, xdiag::OpSum &H);

int main()
try
{
    // Read OpSum from file
    xdiag::OpSum ops = read_opsum(fl, "Interactions");

    // Defines exchange coupling strenght -- This ratio of couplings corresponds to the dimmer phase.
    ops["Jd"] = 1.0;
    ops["J"] = 0.630;

    // Linear Array with target temperatures
    arma::vec Temp = arma::linspace<arma::vec>(0.01, 0.35, 64);

    // Array to store the specific heat
    arma::mat Obser = arma::mat(Temp.n_elem, 3, arma::fill::zeros);

    Tqp_mag_sector(Obser, Temp, ops);

    auto save_fl = FileH5("shastry_sutherland_L_5_W_4.h5", "w!");
    save_fl["Temp"] = Temp;
    save_fl["Observable"] = Obser;

    return 0;
}
catch (Error e)
{
    error_trace(e);
}

void Tqp_mag_sector(arma::mat &Obser, arma::vec &Temp, xdiag::OpSum &H)
{
    auto block = Spinhalf(Nsites); // Create spin-1/2 block with conservation of Sz

    arma::vec eigs;
    arma::mat vecs;
    for (unsigned k = 0; k < Rtqp; k++) // #Perform the calculation for each random vector
    {
        auto res = eigvals_lanczos(H, block, 1, 1e-12, 150,1e-7,k); // Perform the Lanczos interation starting from a random vector;
        
        //This part can be done in the post-processing        
        arma::mat tmatrix = arma::diagmat(res.alphas);
        tmatrix += arma::diagmat(res.betas.head(res.betas.size() - 1), 1) + arma::diagmat(res.betas.head(res.betas.size() - 1), -1);
        
        arma::eig_sym(eigs, vecs, tmatrix); // get eigenvector of T

        for (unsigned k = 0; k < Temp.n_elem; k++)
        {
            arma::vec psi1 =(arma::exp(-(eigs-eigs(0)) / (2.0*Temp(k))) % arma::conj((vecs.row(0).t()))); // compute psi_1

            Obser(k, 0) += arma::dot(psi1, psi1); // get norm factor for the partition function

            Obser(k, 1) += arma::dot(psi1, eigs % psi1); // measure energy

            Obser(k, 2) += arma::dot(psi1,eigs % eigs% psi1); // measure energy square
        }
    }
}