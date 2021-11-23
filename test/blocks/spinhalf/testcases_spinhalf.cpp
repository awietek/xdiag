#include "testcases_spinhalf.h"

#include <cmath>

namespace hydra::testcases::spinhalf {

std::tuple<BondList, Couplings> HBchain(int n_sites, double J1, double J2) {

  BondList bondlist;
  Couplings couplings;
  couplings["J1"] = J1;
  couplings["J2"] = J2;

  for (int s = 0; s < n_sites; ++s) {
    bondlist << Bond("HB", "J1", {s, (s + 1) % n_sites});
  }
  for (int s = 0; s < n_sites; ++s) {
    bondlist << Bond("HB", "J2", {s, (s + 2) % n_sites});
  }
  return {bondlist, couplings};
}

std::tuple<BondList, Couplings, lila::Vector<double>>
HBchain_fullspectrum_nup(int L, int nup) {

  auto [bondlist, couplings] = HBchain(L, 1.0);

  lila::Vector<double> eigs;
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
  return {bondlist, couplings, eigs};
}

std::tuple<BondList, Couplings> HB_alltoall(int n_sites) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  BondList bondlist;
  Couplings couplings;
  for (int s1 = 0; s1 < n_sites; ++s1)
    for (int s2 = s1 + 1; s2 < n_sites; ++s2) {
      std::stringstream ss;
      ss << "J" << s1 << "_" << s2;
      std::string name = ss.str();
      double value = distribution(generator);
      bondlist << Bond("HB", name, {s1, s2});
      couplings[name] = value;
    }
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings, double> triangular_12_complex(int nup,
                                                              double eta) {

  BondList bonds;
  bonds << Bond("ISING", "Jz", {0, 5});
  bonds << Bond("ISING", "Jz", {8, 2});
  bonds << Bond("ISING", "Jz", {2, 7});
  bonds << Bond("ISING", "Jz", {1, 4});
  bonds << Bond("ISING", "Jz", {4, 8});
  bonds << Bond("ISING", "Jz", {6, 10});
  bonds << Bond("ISING", "Jz", {10, 0});
  bonds << Bond("ISING", "Jz", {5, 9});
  bonds << Bond("ISING", "Jz", {9, 3});
  bonds << Bond("ISING", "Jz", {3, 6});
  bonds << Bond("ISING", "Jz", {7, 11});
  bonds << Bond("ISING", "Jz", {11, 1});
  bonds << Bond("ISING", "Jz", {0, 8});
  bonds << Bond("ISING", "Jz", {8, 6});
  bonds << Bond("ISING", "Jz", {2, 10});
  bonds << Bond("ISING", "Jz", {1, 9});
  bonds << Bond("ISING", "Jz", {4, 3});
  bonds << Bond("ISING", "Jz", {6, 1});
  bonds << Bond("ISING", "Jz", {10, 4});
  bonds << Bond("ISING", "Jz", {5, 2});
  bonds << Bond("ISING", "Jz", {9, 7});
  bonds << Bond("ISING", "Jz", {3, 11});
  bonds << Bond("ISING", "Jz", {7, 0});
  bonds << Bond("ISING", "Jz", {11, 5});
  bonds << Bond("ISING", "Jz", {0, 4});
  bonds << Bond("ISING", "Jz", {8, 3});
  bonds << Bond("ISING", "Jz", {2, 6});
  bonds << Bond("ISING", "Jz", {1, 5});
  bonds << Bond("ISING", "Jz", {4, 9});
  bonds << Bond("ISING", "Jz", {6, 11});
  bonds << Bond("ISING", "Jz", {10, 1});
  bonds << Bond("ISING", "Jz", {5, 8});
  bonds << Bond("ISING", "Jz", {9, 2});
  bonds << Bond("ISING", "Jz", {3, 7});
  bonds << Bond("ISING", "Jz", {7, 10});
  bonds << Bond("ISING", "Jz", {11, 0});
  bonds << Bond("EXCHANGE", "Jx", {0, 8});
  bonds << Bond("EXCHANGE", "Jx", {8, 6});
  bonds << Bond("EXCHANGE", "Jx", {2, 10});
  bonds << Bond("EXCHANGE", "Jx", {1, 9});
  bonds << Bond("EXCHANGE", "Jx", {4, 3});
  bonds << Bond("EXCHANGE", "Jx", {6, 1});
  bonds << Bond("EXCHANGE", "Jx", {10, 4});
  bonds << Bond("EXCHANGE", "Jx", {5, 2});
  bonds << Bond("EXCHANGE", "Jx", {9, 7});
  bonds << Bond("EXCHANGE", "Jx", {3, 11});
  bonds << Bond("EXCHANGE", "Jx", {7, 0});
  bonds << Bond("EXCHANGE", "Jx", {11, 5});
  bonds << Bond("EXCHANGE", "Jx", {0, 10});
  bonds << Bond("EXCHANGE", "Jx", {8, 4});
  bonds << Bond("EXCHANGE", "Jx", {2, 8});
  bonds << Bond("EXCHANGE", "Jx", {1, 11});
  bonds << Bond("EXCHANGE", "Jx", {4, 1});
  bonds << Bond("EXCHANGE", "Jx", {6, 3});
  bonds << Bond("EXCHANGE", "Jx", {10, 6});
  bonds << Bond("EXCHANGE", "Jx", {5, 0});
  bonds << Bond("EXCHANGE", "Jx", {9, 5});
  bonds << Bond("EXCHANGE", "Jx", {3, 9});
  bonds << Bond("EXCHANGE", "Jx", {7, 2});
  bonds << Bond("EXCHANGE", "Jx", {11, 7});
  bonds << Bond("EXCHANGE", "Jx", {0, 11});
  bonds << Bond("EXCHANGE", "Jx", {8, 5});
  bonds << Bond("EXCHANGE", "Jx", {2, 9});
  bonds << Bond("EXCHANGE", "Jx", {1, 10});
  bonds << Bond("EXCHANGE", "Jx", {4, 0});
  bonds << Bond("EXCHANGE", "Jx", {6, 2});
  bonds << Bond("EXCHANGE", "Jx", {10, 7});
  bonds << Bond("EXCHANGE", "Jx", {5, 1});
  bonds << Bond("EXCHANGE", "Jx", {9, 4});
  bonds << Bond("EXCHANGE", "Jx", {3, 8});
  bonds << Bond("EXCHANGE", "Jx", {7, 3});
  bonds << Bond("EXCHANGE", "Jx", {11, 6});

  Couplings cpls;
  cpls["Jz"] = 1.0;
  cpls["Jx"] = complex(cos(2 * M_PI * eta), sin(2 * M_PI * eta));
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
  return {bonds, cpls, e0};
}

} // namespace hydra::testcases::spinhalf
