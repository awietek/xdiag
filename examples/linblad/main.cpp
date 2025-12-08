#include <xdiag/all.hpp>
using namespace xdiag;
int main()
try
{
    // system parameters:
    int L = 6; // total number of sites
    double gamma = 2.0; // dissipation strength
    int Nsites_liouville_space = 2 * L; // the effective Hilbert Space is doubled!
    int nup = Nsites_liouville_space / 2; // S^z = 0 sector
    
    // build OpSum
    arma::cx_mat sx = arma::cx_mat({{0, std::complex<double>(0.5,0.0)},{std::complex<double>(0.5,0.0), 0}});
    arma::cx_mat sy = arma::cx_mat({{0, std::complex<double>(0.0,-0.5)},{std::complex<double>(0.0,0.5), 0}});

    arma::cx_mat sx_sx = arma::kron(sx,sx);
    arma::cx_mat sy_sy = arma::kron(sy,sy);

    //
    Spinhalf block(Nsites_liouville_space, nup);
    
    OpSum ops;
    for (int i = 0; i < L; ++i) {
        ops += std::complex<double>(0.0,-1.0)* Op("Matrix", {i, (i + 1) % L},sx_sx);
        ops += std::complex<double>(0.0,-1.0)* Op("Matrix", {i, (i + 1) % L},sy_sy);

        ops += std::complex<double>(0.0,1.0)* Op("Matrix", {L+i,L + (i + 1) % L},sx_sx);
        ops += std::complex<double>(0.0,1.0)* Op("Matrix", {L+i,L + (i + 1) % L},sy_sy);
        ops += gamma* Op("SzSz", {i, L + i});
    }

    // convert to matrix to do Full ED:
    arma::cx_vec eigs;

    arma::cx_mat Lmat = matrixC(ops, block); // convert to matrix
    arma::eig_gen(eigs, Lmat);       // perform exact diagonalization

    eigs -= L*gamma/4.0; // correct zero point value

    auto save_fl = FileH5("bulk_dephasing_spectrum.h5", "w!");
    save_fl["eigenvalues"] = eigs;
    
    return 0;
}
catch (Error e)
{
    error_trace(e);
}