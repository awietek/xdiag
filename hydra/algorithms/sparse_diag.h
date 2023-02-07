#pragma once

#include <utility>

#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

namespace hydra {

// Computing the ground state energy
double eig0_real(BondList const &bondlist, Block const &block,
                 double precision = 1e-12, int max_iterations = 1000,
                 uint64_t seed = 42);

double eig0_cplx(BondList const &bondlist, Block const &block,
                 double precision = 1e-12, int max_iterations = 1000,
                 uint64_t seed = 42);

double eig0(BondList const &bondlist, Block const &block,
            double precision = 1e-12, int max_iterations = 1000,
            uint64_t seed = 42);

// Computing the ground state
StateReal groundstate_real(BondList const &bondlist, Block const &block,
                           double precision = 1e-12, int max_iterations = 1000,
                           uint64_t seed = 42);
StateCplx groundstate_cplx(BondList const &bondlist, Block const &block,
                           double precision = 1e-12, int max_iterations = 1000,
                           uint64_t seed = 42);

StateCplx groundstate(BondList const &bondlist, Block const &block,
                      double precision = 1e-12, int max_iterations = 1000,
                      uint64_t seed = 42);

// Computing both the ground state and its energy
std::pair<double, StateReal> eig0_groundstate_real(BondList const &bondlist,
                                                   Block const &block,
                                                   double precision = 1e-12,
                                                   int max_iterations = 1000,
                                                   uint64_t seed = 42);

std::pair<double, StateCplx> eig0_groundstate_cplx(BondList const &bondlist,
                                                   Block const &block,
                                                   double precision = 1e-12,
                                                   int max_iterations = 1000,
                                                   uint64_t seed = 42);

std::pair<double, StateCplx> eig0_groundstate(BondList const &bondlist,
                                              Block const &block,
                                              double precision = 1e-12,
                                              int max_iterations = 1000,
                                              uint64_t seed = 42);

} // namespace hydra
