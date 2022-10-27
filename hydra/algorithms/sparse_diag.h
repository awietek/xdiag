#pragma once

#include <utility>

#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

// Computing the ground state energy
double eig0_real(BondList const &bondlist, Block const &block,
                 double precision = 1e-12, int seed = 42,
                 int max_iterations = 1000);

double eig0_cplx(BondList const &bondlist, Block const &block,
                 double precision = 1e-12, int seed = 42,
                 int max_iterations = 1000);

double eig0(BondList const &bondlist, Block const &block,
            double precision = 1e-12, int seed = 42, int max_iterations = 1000);

// Computing the ground state
StateReal groundstate_real(BondList const &bondlist, Block const &block,
                           double precision = 1e-12, int seed = 42,
                           int max_iterations = 1000);
StateCplx groundstate_cplx(BondList const &bondlist, Block const &block,
                           double precision = 1e-12, int seed = 42,
                           int max_iterations = 1000);

StateCplx groundstate(BondList const &bondlist, Block const &block,
                      double precision = 1e-12, int seed = 42,
                      int max_iterations = 1000);

// Computing both the ground state and its energy
std::pair<double, StateReal> eig0_groundstate_real(BondList const &bondlist,
                                                   Block const &block,
                                                   double precision = 1e-12,
                                                   int seed = 42,
                                                   int max_iterations = 1000);

std::pair<double, StateCplx> eig0_groundstate_cplx(BondList const &bondlist,
                                                   Block const &block,
                                                   double precision = 1e-12,
                                                   int seed = 42,
                                                   int max_iterations = 1000);

std::pair<double, StateCplx> eig0_groundstate(BondList const &bondlist,
                                              Block const &block,
                                              double precision = 1e-12,
                                              int seed = 42,
                                              int max_iterations = 1000);

} // namespace hydra
