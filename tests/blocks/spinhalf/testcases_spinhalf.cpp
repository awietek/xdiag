#include "testcases_spinhalf.hpp"

#include <cmath>

namespace xdiag::testcases::spinhalf {

OpSum HBchain(int64_t n_sites, double J1, double J2) {

  OpSum ops;
  ops["J1"] = J1;
  ops["J2"] = J2;

  for (int64_t s = 0; s < n_sites; ++s) {
    ops += "J1" * Op("SdotS", {s, (s + 1) % n_sites});
  }
  if (n_sites > 2) {
    for (int64_t s = 0; s < n_sites; ++s) {
      ops += "J2" * Op("SdotS", {s, (s + 2) % n_sites});
    }
  }
  return ops;
}

std::tuple<OpSum, arma::Col<double>> HBchain_fullspectrum_nup(int64_t L,
                                                              int64_t nup) {

  auto ops = HBchain(L, 1.0);

  arma::Col<double> eigs;
  if (L == 2) {
    if ((nup == 0) || (nup == 2))
      eigs = {0.5};
    else if (nup == 1)
      eigs = {-1.5, 0.5};

  } else if (L == 3) {
    if ((nup == 0) || (nup == 3))
      eigs = {0.75};
    else if ((nup == 1) || (nup == 2))
      eigs = {-0.75, -0.75, 0.75};

  } else if (L == 4) {
    if (nup == 2)
      eigs = {-2.0, -1.0, 0.0, 0.0, 0.0, 1.0};
    else if ((nup == 1) || (nup == 3))
      eigs = {-1.0, 0.0, 0.0, 1.0};
    else if ((nup == 0) || (nup == 4))
      eigs = {1.0};

  } else if (L == 5) {
    if ((nup == 2) || (nup == 3)) {
      eigs = {-1.868033988749894,
              -1.868033988749894,
              -0.75,
              -0.559016994374947,
              -0.559016994374947,
              0.368033988749894,
              0.368033988749895,
              0.559016994374947,
              0.559016994374947,
              1.25};
    } else if ((nup == 1) || (nup == 4)) {
      eigs = {-0.559016994374948, -0.559016994374947, 0.559016994374947,
              0.559016994374947, 1.25};
    } else if ((nup == 0) || (nup == 5)) {
      eigs = {1.25};
    }

  } else if (L == 6) {
    if (nup == 3) {
      eigs = {-2.802775637731995e+00, -2.118033988749895e+00,
              -1.499999999999999e+00, -1.280776406404416e+00,
              -1.280776406404415e+00, -1.000000000000000e+00,
              -1.000000000000000e+00, -5.000000000000003e-01,
              -5.000000000000002e-01, -5.000000000000000e-01,
              -1.082836864502427e-16, 1.333450508964792e-16,
              1.180339887498945e-01,  4.999999999999998e-01,
              7.807764064044148e-01,  7.807764064044157e-01,
              8.027756377319950e-01,  1.000000000000001e+00,
              1.000000000000001e+00,  1.499999999999999e+00};
    } else if ((nup == 2) || (nup == 4)) {
      eigs = {-2.118033988749896e+00, -1.280776406404416e+00,
              -1.280776406404415e+00, -1.000000000000000e+00,
              -9.999999999999999e-01, -5.000000000000000e-01,
              -1.602135480420095e-16, -1.020727320943357e-16,
              1.180339887498948e-01,  4.999999999999999e-01,
              7.807764064044147e-01,  7.807764064044155e-01,
              9.999999999999997e-01,  1.000000000000000e+00,
              1.500000000000000e+00};
    } else if ((nup == 1) || (nup == 5)) {
      eigs = {-5.000000000000004e-01, 2.255027837597781e-17,
              3.556194393869664e-16,  9.999999999999998e-01,
              9.999999999999999e-01,  1.500000000000000e+00};
    } else if ((nup == 0) || (nup == 6)) {
      eigs = {1.5};
    }
  }
  return {ops, eigs};
}

OpSum HB_alltoall(int64_t n_sites) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  OpSum ops;
  for (int64_t s1 = 0; s1 < n_sites; ++s1)
    for (int64_t s2 = s1 + 1; s2 < n_sites; ++s2) {
      std::stringstream ss;
      ss << "J" << s1 << "_" << s2;
      std::string name = ss.str();
      double value = distribution(generator);
      ops += name * Op("SdotS", {s1, s2});
      ops[name] = value;
    }
  return ops;
}

std::tuple<OpSum, double> triangular_12_complex(int64_t nup, double eta) {

  OpSum ops;
  ops += "Jz" * Op("SzSz", {0, 5});
  ops += "Jz" * Op("SzSz", {8, 2});
  ops += "Jz" * Op("SzSz", {2, 7});
  ops += "Jz" * Op("SzSz", {1, 4});
  ops += "Jz" * Op("SzSz", {4, 8});
  ops += "Jz" * Op("SzSz", {6, 10});
  ops += "Jz" * Op("SzSz", {10, 0});
  ops += "Jz" * Op("SzSz", {5, 9});
  ops += "Jz" * Op("SzSz", {9, 3});
  ops += "Jz" * Op("SzSz", {3, 6});
  ops += "Jz" * Op("SzSz", {7, 11});
  ops += "Jz" * Op("SzSz", {11, 1});
  ops += "Jz" * Op("SzSz", {0, 8});
  ops += "Jz" * Op("SzSz", {8, 6});
  ops += "Jz" * Op("SzSz", {2, 10});
  ops += "Jz" * Op("SzSz", {1, 9});
  ops += "Jz" * Op("SzSz", {4, 3});
  ops += "Jz" * Op("SzSz", {6, 1});
  ops += "Jz" * Op("SzSz", {10, 4});
  ops += "Jz" * Op("SzSz", {5, 2});
  ops += "Jz" * Op("SzSz", {9, 7});
  ops += "Jz" * Op("SzSz", {3, 11});
  ops += "Jz" * Op("SzSz", {7, 0});
  ops += "Jz" * Op("SzSz", {11, 5});
  ops += "Jz" * Op("SzSz", {0, 4});
  ops += "Jz" * Op("SzSz", {8, 3});
  ops += "Jz" * Op("SzSz", {2, 6});
  ops += "Jz" * Op("SzSz", {1, 5});
  ops += "Jz" * Op("SzSz", {4, 9});
  ops += "Jz" * Op("SzSz", {6, 11});
  ops += "Jz" * Op("SzSz", {10, 1});
  ops += "Jz" * Op("SzSz", {5, 8});
  ops += "Jz" * Op("SzSz", {9, 2});
  ops += "Jz" * Op("SzSz", {3, 7});
  ops += "Jz" * Op("SzSz", {7, 10});
  ops += "Jz" * Op("SzSz", {11, 0});
  ops += "Jx" * Op("Exchange", {0, 8});
  ops += "Jx" * Op("Exchange", {8, 6});
  ops += "Jx" * Op("Exchange", {2, 10});
  ops += "Jx" * Op("Exchange", {1, 9});
  ops += "Jx" * Op("Exchange", {4, 3});
  ops += "Jx" * Op("Exchange", {6, 1});
  ops += "Jx" * Op("Exchange", {10, 4});
  ops += "Jx" * Op("Exchange", {5, 2});
  ops += "Jx" * Op("Exchange", {9, 7});
  ops += "Jx" * Op("Exchange", {3, 11});
  ops += "Jx" * Op("Exchange", {7, 0});
  ops += "Jx" * Op("Exchange", {11, 5});
  ops += "Jx" * Op("Exchange", {0, 10});
  ops += "Jx" * Op("Exchange", {8, 4});
  ops += "Jx" * Op("Exchange", {2, 8});
  ops += "Jx" * Op("Exchange", {1, 11});
  ops += "Jx" * Op("Exchange", {4, 1});
  ops += "Jx" * Op("Exchange", {6, 3});
  ops += "Jx" * Op("Exchange", {10, 6});
  ops += "Jx" * Op("Exchange", {5, 0});
  ops += "Jx" * Op("Exchange", {9, 5});
  ops += "Jx" * Op("Exchange", {3, 9});
  ops += "Jx" * Op("Exchange", {7, 2});
  ops += "Jx" * Op("Exchange", {11, 7});
  ops += "Jx" * Op("Exchange", {0, 11});
  ops += "Jx" * Op("Exchange", {8, 5});
  ops += "Jx" * Op("Exchange", {2, 9});
  ops += "Jx" * Op("Exchange", {1, 10});
  ops += "Jx" * Op("Exchange", {4, 0});
  ops += "Jx" * Op("Exchange", {6, 2});
  ops += "Jx" * Op("Exchange", {10, 7});
  ops += "Jx" * Op("Exchange", {5, 1});
  ops += "Jx" * Op("Exchange", {9, 4});
  ops += "Jx" * Op("Exchange", {3, 8});
  ops += "Jx" * Op("Exchange", {7, 3});
  ops += "Jx" * Op("Exchange", {11, 6});

  ops["Jz"] = 1.0;
  ops["Jx"] = complex(cos(2 * M_PI * eta), sin(2 * M_PI * eta));
  // lila::Log("Jx {} {}", lila::real(cpls["Jx"]), lila::imag(cpls["Jx"]));

  double e0 = 0;
  if (nup == 6) {
    if (eta == 0.00) {
      e0 = -7.3239616252546069219;
    } else if (eta == 0.01) {
      e0 = -7.4775439738976006154;
    } else if (eta == 0.02) {
      e0 = -7.810588969572918927;
    } else if (eta == 0.03) {
      e0 = -8.1986453574035831338;
    } else if (eta == 0.04) {
      e0 = -8.5966264674625048059;
    } else if (eta == 0.05) {
      e0 = -8.9863396883370398882;
    }
  }
  return {ops, e0};
}

} // namespace xdiag::testcases::spinhalf
