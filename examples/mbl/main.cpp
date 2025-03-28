#include <xdiag/all.hpp>

using namespace xdiag;

double level_spacing_ratio(arma::vec eigs)
{
    arma::vec diff=arma::vec(eigs.size()-1);
    for (unsigned i=0; i<eigs.size()-1; i++)
    {
        diff[i]=eigs[i+1]-eigs[i];
    }
    arma::vec r = arma::vec(diff.size()-1);
    for (unsigned i=0; i<diff.size()-1; i++)
    {
        r[i]=std::min(diff[i],diff[i+1])/std::max(diff[i],diff[i+1]);
    }
    return arma::sum(r)/r.size();
}


double ipr_func(arma::vec state)
{
    return arma::sum(arma::pow(arma::abs(state),4));
}


int main()
{
    unsigned L = 10;
    unsigned nres = 2;
    unsigned neigs = 4;
    arma::vec Ws = arma::linspace(0,5,100);
    arma::mat rs = arma::mat(nres,Ws.size());
    arma::cube iprs = arma::cube(nres,neigs,Ws.size());

    for (unsigned i=0; i<nres; i++)
    {
        arma::vec J = arma::randu(L-1,arma::distr_param(0.8,1.2));
        for (unsigned j=0; j<Ws.size(); j++)
        {
            arma::vec h = arma::randu(L,arma::distr_param(-Ws(j),Ws(j)));
            // construct Hamiltonian and Hilbert space
            auto block = Spinhalf(L);
            auto ops = OpSum();
            for (unsigned l=0; l<L-1; l++)
            {
                ops += J(l) * Op("SzSz",{l,l+1});
            }
            for (unsigned l=0; l<L; l++)
            {
                ops += h(l) * Op("Sz",l);
                ops += 0.5 * Op("S+",l);
                ops += 0.5 * Op("S-",l);
            }

            arma::mat H = matrix(ops, block);
            arma::mat eigvecs;
            arma::vec eigvals;
            arma::eig_sym(eigvals,eigvecs,H);

            rs(i,j) = level_spacing_ratio(eigvals);

            for (unsigned k=0; k<neigs; k++)
            {
                arma::vec eigvec = eigvecs.col(k);
                iprs(i,j,k) = ipr_func(eigvec);
            }

        }
    }

    arma::vec r = arma::sum(rs,0);
    arma::mat ipr = arma::sum(iprs,0);

    // save and plot

    return 0;
}