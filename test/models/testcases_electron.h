#pragma once

#include <hydra/all.h>
#include <random>

namespace hydra::testcases::electron {

std::tuple<BondList, Couplings> get_linear_chain(int n_sites, double t,
                                                 double U);

template <class bit_t>
std::tuple<SpaceGroup<bit_t>, std::vector<Representation>>
get_cyclic_group_irreps(int n_sites);

template <class bit_t>
std::tuple<SpaceGroup<bit_t>, std::vector<Representation>, std::vector<int>>
get_cyclic_group_irreps_mult(int n_sites);

std::tuple<BondList, Couplings> heisenberg_triangle();

std::tuple<BondList, Couplings> heisenberg_alltoall(int n_sites);

std::tuple<BondList, Couplings> heisenberg_kagome15();
std::tuple<BondList, Couplings> heisenberg_kagome39();

std::tuple<BondList, Couplings> freefermion_alltoall(int n_sites);

std::tuple<BondList, Couplings> freefermion_alltoall_complex_updn(int n_sites);

std::tuple<BondList, Couplings> tJchain(int n_sites, double t, double J);

std::tuple<BondList, Couplings, lila::Vector<double>> randomAlltoAll4NoU();
std::tuple<BondList, Couplings, lila::Vector<double>> randomAlltoAll4();

std::tuple<BondList, Couplings> randomAlltoAll3();

std::tuple<BondList, Couplings> square2x2(double t, double J);

std::tuple<BondList, Couplings> square3x3(double t, double J);

} // namespace hydra::testcases::electron
