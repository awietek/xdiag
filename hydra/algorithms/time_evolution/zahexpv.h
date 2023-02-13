#pragma once
//
// Created by Luke Staszewski on 27.01.23.
//

#include <extern/armadillo/armadillo>
#include <functional>
#include <tuple>

namespace hydra {

std::tuple<double, double>
zahexpv(double time,
        std::function<arma::cx_vec(arma::cx_vec const &)> const &apply_A,
        arma::cx_vec &v, double tol = 1e-12, int m = 5, double anorm = 0.,
        int nnorm = 2);

}
