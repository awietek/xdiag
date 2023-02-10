//
// Created by Luke Staszewski on 07.02.23.
//
#pragma once
#include <cmath>
#include <algorithm>
#include <armadillo>

template <class coeff_t>
inline double log2abs(coeff_t x) {
    double value;

    if (x == 0.0)
        value = -1e30;
    else
        value = std::log(std::abs(x)) / std::log(2.0);

    return value;
}

inline arma::Mat<arma::cx_double> ExpM(arma::Mat<arma::cx_double> const &A, arma::cx_double alpha = 1.) {
    int n = A.n_rows;
    assert(n == A.n_cols);

    const int q = 6;
    arma::Mat<arma::cx_double> A2 = alpha * A;
    auto a_norm = arma::norm(A2, "inf");
    int ee = (int)log2abs(a_norm) + 1;
    int s = std::max(0, ee + 1);
    double t = 1.0 / pow(2.0, s);

    A2 = t*A2;
    arma::Mat<arma::cx_double> X = A2;
    double c = 0.5;
    arma::Mat<arma::cx_double> E = arma::conv_to<arma::Mat<arma::cx_double>>::from(arma::eye(n, n));
    E += c * A2;
    arma::Mat<arma::cx_double> D = arma::conv_to<arma::Mat<arma::cx_double>>::from(arma::eye(n, n));
    D -= c * A2;

    for (int k = 2; k <= q; k++) {
        c *= (q - k + 1) / (double)(k * (2 * q - k + 1));
        X = X*A2;
        E += c * X;

        if (k%2 == 0)
            D += c * X;
        else
            D -= c * X;
    }

    // inverse of D
    E = E * D.i();
    E = arma::powmat(E, std::pow(2, s));

    return E;
}


