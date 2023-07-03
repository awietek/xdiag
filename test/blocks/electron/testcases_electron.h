#pragma once

#include <random>
#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::testcases::electron {

BondList get_linear_chain(int n_sites, double t, double U);

BondList get_linear_chain_hb(int n_sites, double J);

std::tuple<PermutationGroup, std::vector<Representation>>
get_cyclic_group_irreps(int n_sites);

std::tuple<PermutationGroup, std::vector<Representation>, std::vector<int>>
get_cyclic_group_irreps_mult(int n_sites);

BondList heisenberg_triangle();

BondList heisenberg_alltoall(int n_sites);

BondList heisenberg_kagome15();
BondList heisenberg_kagome39();

BondList freefermion_alltoall(int n_sites);

BondList freefermion_alltoall_complex_updn(int n_sites);

// BondList tJchain(int n_sites, double t, double J);

std::tuple<BondList, arma::Col<double>> randomAlltoAll4NoU();
std::tuple<BondList, arma::Col<double>> randomAlltoAll4();

BondList randomAlltoAll3();

BondList square2x2(double t, double J);

BondList square3x3(double t, double J);

} // namespace hydra::testcases::electron
