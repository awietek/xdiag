//
// Created by Luke Staszewski on 07.02.23.
//

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include "../pade_matrix_exponential.hpp"
#include <complex>
using namespace std;
using namespace arma;


TEST_CASE("exp of identity", "[pade]"){

    auto I = eye(3,3);
    cx_double i = {1, 0};
    Mat<cx_double> B = i*I;
    B.print();
    Mat<cx_double> E = ExpM(B);
    E.print();

}