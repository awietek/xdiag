//
// Created by Luke Staszewski on 08.02.23.
//
#include "lanczos_time_evolve.hpp"
#include <chrono>
#include "armadillo"

using namespace std;
using namespace arma;
using namespace xdiag;

int main(){
    cout << "timing the lanczos algorithm ... ";

    // defining a model with xdiag
    int n_sites, nup, ndn;
    n_sites = 10; nup = 3; ndn = 3;
    double t, t1, U;
    t = 1; t1 = .5; U = 4;

    auto block = Electron(n_sites, nup, ndn);
    BondList bonds;
    for(int i = 0; i< n_sites; i++){
        bonds << Bond("HOP", "t", {i, (i+1)%n_sites});
        bonds << Bond("HOP", "t1", {i, (i+2)%n_sites});
    }
    bonds["t"] = t; bonds["t1"] = t1; bonds["U"] = U;

    auto psi_0_list = ProductState({"Up", "Up", "Emp", "UpDn", "Emp", "Dn", "Dn", "Emp" , "Emp", "Emp"});
    auto psi_0 =  State(block, psi_0_list);

    // doing the time evolution
    double time = 30;
    auto tic = chrono::high_resolution_clock::now();
    auto psi = zahexpv(time, bonds, psi_0, 1e-2, 10);
    auto toc = chrono::high_resolution_clock::now();
    cout << norm(psi.vector()) << endl;
    cout << "time taken: " <<  (toc - tic).count()*1e-9 << endl;

    // timing with last xdiag

    auto tic1 = chrono::high_resolution_clock::now();
    auto psi_1 = time_evolve(bonds, psi_0, time, 1e-2);
    auto toc1 = chrono::high_resolution_clock::now();
    cout << norm(psi.vector()) << endl;
    cout << "time taken xdiag: " <<  (toc1 - tic1).count()*1e-9 << endl;
    return 0;
}





