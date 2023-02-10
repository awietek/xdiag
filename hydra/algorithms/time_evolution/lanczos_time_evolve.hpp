//
// Created by Luke Staszewski on 27.01.23.
//
#pragma once
#ifndef LANC_TIME_EVOLUTION_TIME_EVOLVE_H
#define LANC_TIME_EVOLUTION_TIME_EVOLVE_H

#endif //LANC_TIME_EVOLUTION_TIME_EVOLVE_H
#include <functional>
#include <armadillo>
#include <hydra/all.h>
using namespace hydra;

arma::Col<arma::cx_double> zahexpv(double time,
                                  const std::function<arma::Col<arma::cx_double> (arma::Col<arma::cx_double> const &)> &apply_A,
                                  arma::Col<arma::cx_double> &v,
                                  double anorm,
                                  double tol,
                                  int m);


template <typename StateT>
StateT zahexpv(double time, BondList const &bonds, StateT &state,
              double tol = 1e-7, int m = 5) {

    auto const &block = state.block();

    arma::cx_double i = {0,1};

    // defining apply_A
    auto mult = [&i, &bonds, &block](arma::Col<arma::cx_double> const &v){
        arma::Col<arma::cx_double> w(v.size(), arma::fill::zeros);
        apply(bonds, block, v, block, w);
        w *= -i;
        return w;
    };

    // getting a norm for the hamiltonian
    arma::Mat<arma::cx_double> H = matrix_cplx(bonds, block);
    double anorm = arma::norm(H, "inf");

    arma::cx_mat H_conj = arma::conj(H);
    auto &v0 = state.vector();

    // actually performing the expokit routine for matrix exponentiation
    arma::Col<arma::cx_double> w = zahexpv(time, mult, v0, anorm, tol, m);

    return State(block, w);
}
