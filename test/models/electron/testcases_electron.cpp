#include "testcases_electron.h"

namespace hydra::testcases::electron {

std::tuple<BondList, Couplings> get_linear_chain(int n_sites, double t,
                                                 double U) {
  // Create model
  BondList bondlist;
  for (int s = 0; s < n_sites; ++s)
    bondlist << Bond("HOP", "T", {s, (s + 1) % n_sites});
  Couplings couplings;
  couplings["T"] = t;
  couplings["U"] = U;
  return {bondlist, couplings};
}

std::tuple<BondList, Couplings> get_linear_chain_hb(int n_sites, double t,
                                                    double U, double J) {
  // Create model
  BondList bondlist;
  for (int s = 0; s < n_sites; ++s)
    bondlist << Bond("HOP", "T", {s, (s + 1) % n_sites});
  for (int s = 0; s < n_sites; ++s)
    bondlist << Bond("HB", "J", {s, (s + 1) % n_sites});

  Couplings couplings;
  couplings["T"] = t;
  couplings["U"] = U;
  couplings["J"] = J;
  return {bondlist, couplings};
}

template <class bit_t>
std::tuple<PermutationGroup, std::vector<Representation>>
get_cyclic_group_irreps(int n_sites) {
  // Create cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = PermutationGroup(permutations);

  // Create irreducible representations
  std::vector<Representation> irreps;
  for (int k = 0; k < n_sites; ++k) {
    std::vector<complex> chis;
    for (int l = 0; l < n_sites; ++l)
      chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                      std::sin(2 * M_PI * l * k / n_sites)});
    auto irrep = Representation(chis);
    irreps.push_back(irrep);
  }
  return {space_group, irreps};
}

template std::tuple<PermutationGroup, std::vector<Representation>>
get_cyclic_group_irreps<uint16>(int n_sites);
template std::tuple<PermutationGroup, std::vector<Representation>>
get_cyclic_group_irreps<uint32>(int n_sites);
template std::tuple<PermutationGroup, std::vector<Representation>>
get_cyclic_group_irreps<uint64>(int n_sites);

template <class bit_t>
std::tuple<PermutationGroup, std::vector<Representation>, std::vector<int>>
get_cyclic_group_irreps_mult(int n_sites) {
  // Create cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = PermutationGroup(permutations);

  // Create irreducible representations
  std::vector<Representation> irreps;
  std::vector<int> multiplicities;
  for (int k = 0; k < n_sites; ++k) {
    std::vector<complex> chis;
    for (int l = 0; l < n_sites; ++l)
      chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                      std::sin(2 * M_PI * l * k / n_sites)});
    auto irrep = Representation(chis);
    irreps.push_back(irrep);
    multiplicities.push_back(1);
  }
  return {space_group, irreps, multiplicities};
}

template std::tuple<PermutationGroup, std::vector<Representation>,
                    std::vector<int>>
get_cyclic_group_irreps_mult<uint16>(int n_sites);
template std::tuple<PermutationGroup, std::vector<Representation>,
                    std::vector<int>>
get_cyclic_group_irreps_mult<uint32>(int n_sites);
template std::tuple<PermutationGroup, std::vector<Representation>,
                    std::vector<int>>
get_cyclic_group_irreps_mult<uint64>(int n_sites);

std::tuple<BondList, Couplings> heisenberg_triangle() {
  BondList bondlist;
  bondlist << Bond("HEISENBERG", "J", {0, 1});
  bondlist << Bond("HEISENBERG", "J", {1, 2});
  bondlist << Bond("HEISENBERG", "J", {2, 0});

  Couplings couplings;
  couplings["J"] = 1.0;
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> heisenberg_alltoall(int n_sites) {
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
      bondlist << Bond("HEISENBERG", name, {s1, s2});
      couplings[name] = value;
    }
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> heisenberg_kagome15() {

  BondList bondlist;
  bondlist << Bond("HEISENBERG", "J", {0, 1});
  bondlist << Bond("HEISENBERG", "J", {0, 5});
  bondlist << Bond("HEISENBERG", "J", {1, 2});
  bondlist << Bond("HEISENBERG", "J", {1, 6});
  bondlist << Bond("HEISENBERG", "J", {2, 3});
  bondlist << Bond("HEISENBERG", "J", {2, 6});
  bondlist << Bond("HEISENBERG", "J", {2, 10});
  bondlist << Bond("HEISENBERG", "J", {3, 4});
  bondlist << Bond("HEISENBERG", "J", {3, 10});
  bondlist << Bond("HEISENBERG", "J", {3, 14});
  bondlist << Bond("HEISENBERG", "J", {4, 5});
  bondlist << Bond("HEISENBERG", "J", {4, 14});
  bondlist << Bond("HEISENBERG", "J", {6, 7});
  bondlist << Bond("HEISENBERG", "J", {7, 8});
  bondlist << Bond("HEISENBERG", "J", {8, 9});
  bondlist << Bond("HEISENBERG", "J", {9, 10});
  bondlist << Bond("HEISENBERG", "J", {9, 11});
  bondlist << Bond("HEISENBERG", "J", {10, 11});
  bondlist << Bond("HEISENBERG", "J", {11, 12});
  bondlist << Bond("HEISENBERG", "J", {12, 13});
  bondlist << Bond("HEISENBERG", "J", {13, 14});
  Couplings couplings;
  couplings["J"] = 1.0;
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> heisenberg_kagome39() {
  BondList bondlist;
  bondlist << Bond("HEISENBERG", "J", {0, 1});
  bondlist << Bond("HEISENBERG", "J", {0, 5});
  bondlist << Bond("HEISENBERG", "J", {0, 27});
  bondlist << Bond("HEISENBERG", "J", {0, 28});
  bondlist << Bond("HEISENBERG", "J", {1, 2});
  bondlist << Bond("HEISENBERG", "J", {1, 6});
  bondlist << Bond("HEISENBERG", "J", {1, 28});
  bondlist << Bond("HEISENBERG", "J", {2, 3});
  bondlist << Bond("HEISENBERG", "J", {2, 6});
  bondlist << Bond("HEISENBERG", "J", {2, 10});
  bondlist << Bond("HEISENBERG", "J", {3, 4});
  bondlist << Bond("HEISENBERG", "J", {3, 10});
  bondlist << Bond("HEISENBERG", "J", {3, 14});
  bondlist << Bond("HEISENBERG", "J", {4, 5});
  bondlist << Bond("HEISENBERG", "J", {4, 14});
  bondlist << Bond("HEISENBERG", "J", {4, 38});
  bondlist << Bond("HEISENBERG", "J", {5, 27});
  bondlist << Bond("HEISENBERG", "J", {5, 38});
  bondlist << Bond("HEISENBERG", "J", {6, 7});
  bondlist << Bond("HEISENBERG", "J", {6, 29});
  bondlist << Bond("HEISENBERG", "J", {7, 8});
  bondlist << Bond("HEISENBERG", "J", {7, 19});
  bondlist << Bond("HEISENBERG", "J", {7, 29});
  bondlist << Bond("HEISENBERG", "J", {8, 9});
  bondlist << Bond("HEISENBERG", "J", {8, 15});
  bondlist << Bond("HEISENBERG", "J", {8, 19});
  bondlist << Bond("HEISENBERG", "J", {9, 10});
  bondlist << Bond("HEISENBERG", "J", {9, 11});
  bondlist << Bond("HEISENBERG", "J", {9, 15});
  bondlist << Bond("HEISENBERG", "J", {10, 11});
  bondlist << Bond("HEISENBERG", "J", {11, 12});
  bondlist << Bond("HEISENBERG", "J", {11, 18});
  bondlist << Bond("HEISENBERG", "J", {12, 13});
  bondlist << Bond("HEISENBERG", "J", {12, 18});
  bondlist << Bond("HEISENBERG", "J", {12, 26});
  bondlist << Bond("HEISENBERG", "J", {13, 14});
  bondlist << Bond("HEISENBERG", "J", {13, 26});
  bondlist << Bond("HEISENBERG", "J", {13, 37});
  bondlist << Bond("HEISENBERG", "J", {14, 37});
  bondlist << Bond("HEISENBERG", "J", {15, 16});
  bondlist << Bond("HEISENBERG", "J", {15, 22});
  bondlist << Bond("HEISENBERG", "J", {16, 17});
  bondlist << Bond("HEISENBERG", "J", {16, 22});
  bondlist << Bond("HEISENBERG", "J", {16, 33});
  bondlist << Bond("HEISENBERG", "J", {17, 18});
  bondlist << Bond("HEISENBERG", "J", {17, 23});
  bondlist << Bond("HEISENBERG", "J", {17, 33});
  bondlist << Bond("HEISENBERG", "J", {18, 23});
  bondlist << Bond("HEISENBERG", "J", {19, 20});
  bondlist << Bond("HEISENBERG", "J", {19, 30});
  bondlist << Bond("HEISENBERG", "J", {20, 21});
  bondlist << Bond("HEISENBERG", "J", {20, 30});
  bondlist << Bond("HEISENBERG", "J", {20, 31});
  bondlist << Bond("HEISENBERG", "J", {21, 22});
  bondlist << Bond("HEISENBERG", "J", {21, 31});
  bondlist << Bond("HEISENBERG", "J", {21, 32});
  bondlist << Bond("HEISENBERG", "J", {22, 32});
  bondlist << Bond("HEISENBERG", "J", {23, 24});
  bondlist << Bond("HEISENBERG", "J", {23, 34});
  bondlist << Bond("HEISENBERG", "J", {24, 25});
  bondlist << Bond("HEISENBERG", "J", {24, 34});
  bondlist << Bond("HEISENBERG", "J", {24, 35});
  bondlist << Bond("HEISENBERG", "J", {25, 26});
  bondlist << Bond("HEISENBERG", "J", {25, 35});
  bondlist << Bond("HEISENBERG", "J", {25, 36});
  bondlist << Bond("HEISENBERG", "J", {26, 36});
  Couplings couplings;
  couplings["J"] = 1.0;
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> freefermion_alltoall(int n_sites) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  BondList bondlist;
  Couplings couplings;
  for (int s1 = 0; s1 < n_sites; ++s1)
    for (int s2 = s1 + 1; s2 < n_sites; ++s2) {
      std::stringstream ss;
      ss << "T" << s1 << "_" << s2;
      std::string name = ss.str();
      double value = distribution(generator);
      bondlist << Bond("HOP", name, {s1, s2});
      couplings[name] = value;
    }
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> freefermion_alltoall_complex_updn(int n_sites) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  BondList bondlist;
  Couplings couplings;
  for (int s1 = 0; s1 < n_sites; ++s1)
    for (int s2 = s1 + 1; s2 < n_sites; ++s2) {

      // Hopping on upspins
      std::stringstream ss_up;
      ss_up << "TUP" << s1 << "_" << s2;
      std::string name_up = ss_up.str();
      complex value_up =
          complex(distribution(generator), distribution(generator));
      bondlist << Bond("HOPUP", name_up, {s1, s2});
      couplings[name_up] = value_up;

      // Hopping on dnspins
      std::stringstream ss_dn;
      ss_dn << "TDN" << s1 << "_" << s2;
      std::string name_dn = ss_dn.str();
      complex value_dn =
          complex(distribution(generator), distribution(generator));
      bondlist << Bond("HOPDN", name_dn, {s1, s2});
      couplings[name_dn] = value_dn;
    }
  return std::make_tuple(bondlist, couplings);
}

// std::tuple<BondList, Couplings> tJchain(int n_sites, double t,
//                                                double J) {

//   BondList bondlist;
//   Couplings couplings;
//   couplings["T"] = t;
//   couplings["J"] = J;
//   for (int s = 0; s < n_sites; ++s) {
//     bondlist << Bond("HUBBARDHOP", "T", {s, (s + 1) % n_sites});
//     bondlist << Bond("HEISENBERG", "J", {s, (s + 1) % n_sites});
//   }
//   return std::make_tuple(bondlist, couplings);
// }

std::tuple<BondList, Couplings, lila::Vector<double>> randomAlltoAll4NoU() {
  BondList bondlist;
  Couplings couplings;
  // couplings["T01"] = 3;
  // couplings["J01"] = 1;
  // couplings["T02"] = 3;
  // couplings["J02"] = -3;
  // couplings["T03"] = 3;
  // couplings["J03"] = 5;
  // couplings["T12"] = 4;
  // couplings["J12"] = -5;
  // couplings["T13"] = -1;
  // couplings["J13"] = -1;
  // couplings["T23"] = 2;
  // couplings["J23"] = 1;

  couplings["T01"] = -3;
  couplings["J01"] = 1;
  couplings["T02"] = -3;
  couplings["J02"] = -3;
  couplings["T03"] = -3;
  couplings["J03"] = 5;
  couplings["T12"] = -4;
  couplings["J12"] = -5;
  couplings["T13"] = 1;
  couplings["J13"] = -1;
  couplings["T23"] = -2;
  couplings["J23"] = 1;

  bondlist << Bond("HUBBARDHOP", "T01", {0, 1});
  bondlist << Bond("HUBBARDHOP", "T02", {0, 2});
  bondlist << Bond("HUBBARDHOP", "T03", {0, 3});
  bondlist << Bond("HUBBARDHOP", "T12", {1, 2});
  bondlist << Bond("HUBBARDHOP", "T13", {1, 3});
  bondlist << Bond("HUBBARDHOP", "T23", {2, 3});
  bondlist << Bond("HEISENBERG", "J01", {0, 1});
  bondlist << Bond("HEISENBERG", "J02", {0, 2});
  bondlist << Bond("HEISENBERG", "J03", {0, 3});
  bondlist << Bond("HEISENBERG", "J12", {1, 2});
  bondlist << Bond("HEISENBERG", "J13", {1, 3});
  bondlist << Bond("HEISENBERG", "J23", {2, 3});

  lila::Vector<double> eigs = {-17.035603173216636,
                               -16.054529653295518,
                               -16.054529653295504,
                               -14.839136196281768,
                               -14.839136196281759,
                               -14.479223672075845,
                               -13.947060439818175,
                               -13.681140962579473,
                               -13.681140962579473,
                               -13.681140962579470,
                               -12.146019505147946,
                               -12.146019505147938,
                               -11.123249987689370,
                               -11.083677166546861,
                               -11.083677166546861,
                               -10.361590604796385,
                               -10.141075725997615,
                               -10.141075725997606,
                               -9.879061771701892,
                               -9.879061771701885,
                               -9.879061771701874,
                               -9.720915055042584,
                               -9.720915055042580,
                               -9.300171000114572,
                               -8.898903149068287,
                               -8.898903149068287,
                               -8.898903149068287,
                               -8.587797030969547,
                               -8.574093646826530,
                               -8.574093646826528,
                               -8.574093646826526,
                               -8.567342877760581,
                               -8.556463239828611,
                               -8.556463239828611,
                               -8.156431544113079,
                               -8.156431544113071,
                               -7.595003505113175,
                               -7.595003505113174,
                               -7.595003505113174,
                               -7.595003505113171,
                               -7.428914803058910,
                               -7.428914803058910,
                               -7.406132446925684,
                               -7.406132446925682,
                               -7.298052445959064,
                               -7.298052445959062,
                               -7.298052445959062,
                               -6.776050147544091,
                               -6.776050147544089,
                               -6.597100597834562,
                               -6.382421301285782,
                               -6.382421301285780,
                               -6.382421301285776,
                               -5.914206919262412,
                               -5.914206919262412,
                               -5.914206919262412,
                               -5.914206919262406,
                               -5.898063094032344,
                               -5.697730676986595,
                               -5.652742708313134,
                               -5.652742708313128,
                               -5.382395669397896,
                               -5.382395669397890,
                               -4.827554533410211,
                               -4.827554533410209,
                               -4.827554533410208,
                               -4.565866945456345,
                               -4.392721098506336,
                               -4.392721098506335,
                               -4.392721098506333,
                               -4.386896721326241,
                               -4.386896721326240,
                               -4.386896721326238,
                               -4.287074157175168,
                               -4.269109370889475,
                               -4.269109370889474,
                               -4.083758285516160,
                               -3.784107949888678,
                               -3.784107949888678,
                               -3.230851175883084,
                               -3.230851175883084,
                               -3.166425888361765,
                               -3.166425888361761,
                               -3.060272421221770,
                               -3.060272421221768,
                               -3.060272421221768,
                               -3.060272421221767,
                               -2.846017856191310,
                               -2.846017856191308,
                               -2.826551366644327,
                               -2.822163676323597,
                               -2.822163676323595,
                               -2.373593337341226,
                               -2.304206358771344,
                               -2.291423386597424,
                               -2.291423386597422,
                               -2.291423386597419,
                               -2.258325746389715,
                               -2.100087802223023,
                               -2.100087802223022,
                               -2.100087802223021,
                               -2.002616246412348,
                               -2.002616246412347,
                               -2.002616246412346,
                               -2.002616246412346,
                               -1.653289828765464,
                               -1.653289828765462,
                               -1.653289828765460,
                               -1.537108454167115,
                               -1.537108454167113,
                               -1.468496478890581,
                               -1.184332042222068,
                               -1.184332042222067,
                               -1.183220245290653,
                               -1.183220245290652,
                               -1.183220245290647,
                               -1.158824284368453,
                               -1.158824284368453,
                               -0.797210829513575,
                               -0.797210829513574,
                               -0.753299251580644,
                               -0.500000000000001,
                               -0.500000000000000,
                               -0.500000000000000,
                               -0.500000000000000,
                               -0.499999999999998,
                               -0.370985460263250,
                               -0.370985460263249,
                               -0.281075696274453,
                               -0.230909105391692,
                               -0.230909105391692,
                               -0.230909105391689,
                               0,
                               0,
                               0.226914386262986,
                               0.226914386262986,
                               0.226914386262988,
                               0.339587764690138,
                               0.339587764690138,
                               0.339587764690140,
                               0.339587764690141,
                               0.864151894040242,
                               0.864151894040242,
                               0.977357729518521,
                               0.977357729518522,
                               0.982508508938287,
                               0.982508508938294,
                               1.184332042222068,
                               1.184332042222068,
                               1.286333102260671,
                               1.360519899915624,
                               1.360519899915626,
                               1.831699701973819,
                               1.831699701973819,
                               1.831699701973821,
                               2.168605503366585,
                               2.304759071083118,
                               2.305593972115476,
                               2.305593972115481,
                               2.305593972115481,
                               2.565136835120275,
                               2.565136835120277,
                               2.680716385503151,
                               2.680716385503155,
                               2.680716385503157,
                               2.859450072401542,
                               2.867740829382918,
                               2.867740829382918,
                               2.867740829382920,
                               2.867740829382922,
                               2.919012177817019,
                               2.919012177817021,
                               3.230851175883083,
                               3.230851175883084,
                               3.586647790757984,
                               3.586647790757985,
                               3.866685809727107,
                               3.866685809727108,
                               3.866685809727108,
                               3.962683310049183,
                               3.962683310049187,
                               3.983903797596342,
                               3.983903797596345,
                               3.983903797596353,
                               4.106914761573067,
                               4.258514587211152,
                               4.258514587211155,
                               4.258514587211158,
                               4.279939892091212,
                               4.647129236685327,
                               4.647129236685331,
                               4.730285398111332,
                               5.382395669397893,
                               5.382395669397895,
                               5.557913081969264,
                               5.729878922142601,
                               5.729878922142602,
                               5.729878922142604,
                               5.729878922142607,
                               5.854994021510809,
                               5.854994021510811,
                               6.026195725670756,
                               6.026195725670764,
                               6.112978522336865,
                               6.112978522336867,
                               6.112978522336871,
                               6.298578032819039,
                               6.627388000300686,
                               6.627388000300687,
                               6.638917394627725,
                               6.638917394627728,
                               6.638917394627730,
                               7.106988282706432,
                               7.261271812957728,
                               7.428914803058909,
                               7.428914803058913,
                               7.634891575794040,
                               7.634891575794041,
                               7.634891575794042,
                               7.634891575794042,
                               8.034109956056216,
                               8.034109956056225,
                               8.433501672445885,
                               8.437627423133117,
                               8.437627423133124,
                               8.437627423133126,
                               8.487415286599031,
                               8.740704187459059,
                               8.740704187459061,
                               8.740704187459061,
                               8.758701982332155,
                               9.740946203547077,
                               9.740946203547077,
                               10.075541640416940,
                               10.075541640416946,
                               10.365553083134904,
                               10.365553083134905,
                               10.898695460947337,
                               10.898695460947337,
                               10.898695460947343,
                               11.368060459508595,
                               11.880069395522252,
                               12.081391252276028,
                               12.081391252276036,
                               12.355338794669144,
                               12.355338794669148,
                               12.833107262067776,
                               14.296824370037875,
                               14.296824370037879,
                               14.296824370037887,
                               15.091839736118505,
                               15.091839736118507,
                               15.880746138642490,
                               17.166681362460483,
                               17.166681362460491,
                               18.194539570876405};

  return std::make_tuple(bondlist, couplings, eigs);
}

std::tuple<BondList, Couplings, lila::Vector<double>> randomAlltoAll4() {
  BondList bondlist;
  Couplings couplings;
  // couplings["U"] = 5;
  // couplings["T01"] = 3;
  // couplings["J01"] = -1;
  // couplings["T02"] = -3;
  // couplings["J02"] = -5;
  // couplings["T03"] = 3;
  // couplings["J03"] = -3;
  // couplings["T12"] = -1;
  // couplings["J12"] = 1;
  // couplings["T13"] = -3;
  // couplings["J13"] = 2;
  // couplings["T23"] = 0;
  // couplings["J23"] = -4;

  couplings["U"] = 5;
  couplings["T01"] = -3;
  couplings["J01"] = -1;
  couplings["T02"] = 3;
  couplings["J02"] = -5;
  couplings["T03"] = -3;
  couplings["J03"] = -3;
  couplings["T12"] = 1;
  couplings["J12"] = 1;
  couplings["T13"] = 3;
  couplings["J13"] = 2;
  couplings["T23"] = 0;
  couplings["J23"] = -4;

  bondlist << Bond("HUBBARDHOP", "T01", {0, 1});
  bondlist << Bond("HUBBARDHOP", "T02", {0, 2});
  bondlist << Bond("HUBBARDHOP", "T03", {0, 3});
  bondlist << Bond("HUBBARDHOP", "T12", {1, 2});
  bondlist << Bond("HUBBARDHOP", "T13", {1, 3});
  bondlist << Bond("HUBBARDHOP", "T23", {2, 3});
  bondlist << Bond("HEISENBERG", "J01", {0, 1});
  bondlist << Bond("HEISENBERG", "J02", {0, 2});
  bondlist << Bond("HEISENBERG", "J03", {0, 3});
  bondlist << Bond("HEISENBERG", "J12", {1, 2});
  bondlist << Bond("HEISENBERG", "J13", {1, 3});
  bondlist << Bond("HEISENBERG", "J23", {2, 3});

  lila::Vector<double> eigs = {
      -12.270231830055396, -12.270231830055389, -10.733666336755952,
      -10.069390063366962, -9.060858377591751,  -9.060858377591751,
      -9.060858377591751,  -8.419560252873284,  -8.419560252873282,
      -8.419560252873278,  -6.383158148575644,  -6.383158148575637,
      -6.383158148575632,  -6.352277902330186,  -6.352277902330185,
      -6.273324224596429,  -6.273324224596422,  -6.250906641891413,
      -6.250906641891411,  -6.164309032262214,  -6.164309032262212,
      -6.164309032262212,  -6.164309032262211,  -5.730618769293929,
      -5.448935789669534,  -5.448935789669532,  -5.038951239070341,
      -4.949876862328434,  -4.949876862328423,  -4.532986251596143,
      -4.532986251596141,  -4.532986251596141,  -4.532986251596141,
      -3.353197524407229,  -3.353197524407228,  -3.353197524407226,
      -3.273406176414287,  -3.002245918852136,  -3.002245918852133,
      -2.753141709527037,  -2.753141709527034,  -2.753141709527034,
      -2.753141709527031,  -2.646622091502864,  -2.646622091502863,
      -2.646622091502862,  -2.500000000000006,  -2.500000000000000,
      -2.500000000000000,  -2.500000000000000,  -2.499999999999995,
      -2.002043720414641,  -2.002043720414640,  -2.002043720414639,
      -1.825844696317418,  -1.825844696317417,  -1.587175599207617,
      -1.587175599207614,  -1.332228906443854,  -1.332228906443853,
      -0.953827936085984,  -0.635382900549627,  -0.635382900549625,
      -0.635382900549624,  -0.397581339114120,  -0.397581339114115,
      -0.302660107585638,  -0.302660107585633,  -0.302660107585631,
      -0.080803683543577,  -0.080803683543570,  0,
      0.216457675593863,   0.256166291525251,   0.601566977837033,
      0.601566977837038,   0.601566977837038,   0.601566977837040,
      0.975606313924293,   0.975606313924297,   1.014605271271656,
      1.014605271271658,   1.015859809070357,   1.015859809070358,
      1.015859809070360,   1.020308313587130,   1.114881844698814,
      1.791357286454801,   1.791357286454802,   1.791357286454810,
      1.812876191672553,   2.051032557261542,   2.051032557261543,
      2.054529439890590,   2.054529439890591,   2.054529439890594,
      2.464728271742040,   2.464728271742042,   2.464728271742044,
      2.464728271742044,   2.561461716067067,   2.599451504679192,
      2.710382274715613,   2.710382274715615,   2.710382274715616,
      2.999999999999997,   3.000000000000001,   3.165957899766594,
      3.165957899766600,   3.217491411604103,   3.217491411604104,
      3.217491411604105,   3.264426167093818,   3.264426167093818,
      3.275854965551124,   3.873065426792698,   3.930285431436003,
      3.930285431436005,   4.357654287008264,   4.373227423701834,
      4.373227423701834,   4.373227423701836,   4.744551703509988,
      4.744551703509988,   4.764857447396031,   4.764857447396040,
      4.764857447396043,   4.838082241099029,   4.838082241099031,
      5.078388983561651,   5.078388983561652,   5.095728306021306,
      5.095728306021313,   5.095728306021315,   5.095728306021317,
      5.270280349321774,   5.629364135391933,   5.629364135391936,
      5.732050357363664,   5.732050357363669,   5.732050357363673,
      5.902527336054253,   5.997898395939853,   5.997898395939854,
      5.997898395939856,   6.168989353312808,   6.168989353312815,
      6.168989353312816,   6.168989353312816,   6.251638235870590,
      6.251638235870590,   6.639239164264768,   6.871779959503020,
      6.871779959503024,   6.913606012729136,   7.197663951269839,
      7.197663951269844,   7.241663577600812,   7.241663577600815,
      7.548559413909176,   7.548559413909178,   7.548559413909183,
      7.889853801872584,   8.012439704238972,   8.012439704238977,
      8.048368645785640,   8.048368645785644,   8.195982486905091,
      8.195982486905091,   8.195982486905095,   8.291793376177347,
      8.291793376177351,   8.291793376177356,   8.468003039901994,
      8.884687492644268,   8.929394188779456,   8.929394188779456,
      9.084392860883399,   9.084392860883403,   9.084392860883410,
      9.119424084472175,   9.119424084472177,   9.119424084472177,
      9.119424084472181,   9.374280903298303,   9.374280903298303,
      9.470513926885971,   9.470513926885971,   9.807459688038790,
      9.894356293199492,   10.161917758900962,  10.161917758900971,
      10.343135676951986,  10.647301560880138,  10.647301560880139,
      10.781521078539114,  10.816967757221121,  10.816967757221125,
      10.816967757221128,  10.989949094629180,  10.989949094629189,
      10.989949094629189,  11.783921524648289,  11.862403712063079,
      11.862403712063086,  11.999999999999995,  11.999999999999995,
      12.122579175754746,  12.122579175754746,  12.422108994367830,
      12.422108994367832,  12.660280744665648,  12.660280744665650,
      12.660280744665654,  12.782275159258591,  13.142554967689829,
      13.262004386769918,  13.262004386769929,  13.262004386769933,
      13.345289206597041,  13.345289206597041,  13.920776472179945,
      14.125358959870058,  14.125358959870061,  14.245875071040452,
      14.245875071040452,  14.710043063865781,  14.821095142124749,
      15.455920358765942,  15.455920358765947,  15.455920358765953,
      15.977688392619838,  15.977688392619839,  16.548872176333433,
      16.587175599207608,  16.587175599207615,  16.668859213157941,
      16.668859213157944,  16.859992272946350,  17.289815197741845,
      17.339978797436935,  17.339978797436938,  17.339978797436956,
      18.075989445512761,  18.075989445512761,  18.524216278708529,
      18.776715574088868,  18.776715574088868,  20.000000000000000,
      20.972807969213903,  21.250906641891415,  21.250906641891415,
      22.411220290848199,  22.411220290848210,  22.508215798228996,
      25.052426347353144};

  return std::make_tuple(bondlist, couplings, eigs);
}

std::tuple<BondList, Couplings> randomAlltoAll3() {
  BondList bondlist;
  Couplings couplings;
  couplings["T01"] = 1;
  couplings["J01"] = -2;
  couplings["T02"] = 0;
  couplings["J02"] = -1;
  couplings["T12"] = -5;
  couplings["J12"] = -3;
  bondlist << Bond("HUBBARDHOP", "T01", {0, 1});
  bondlist << Bond("HUBBARDHOP", "T02", {0, 2});
  bondlist << Bond("HUBBARDHOP", "T12", {1, 2});
  bondlist << Bond("HEISENBERG", "J01", {0, 1});
  bondlist << Bond("HEISENBERG", "J02", {0, 2});
  bondlist << Bond("HEISENBERG", "J12", {1, 2});
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> square2x2(double t, double J) {
  BondList bondlist;
  Couplings couplings;
  couplings["T"] = t;
  couplings["J"] = J;
  bondlist << Bond("HUBBARDHOP", "T", {0, 1});
  bondlist << Bond("HUBBARDHOP", "T", {1, 0});
  bondlist << Bond("HUBBARDHOP", "T", {2, 3});
  bondlist << Bond("HUBBARDHOP", "T", {3, 2});
  bondlist << Bond("HUBBARDHOP", "T", {0, 2});
  bondlist << Bond("HUBBARDHOP", "T", {2, 0});
  bondlist << Bond("HUBBARDHOP", "T", {1, 3});
  bondlist << Bond("HUBBARDHOP", "T", {3, 1});
  bondlist << Bond("HEISENBERG", "J", {0, 1});
  bondlist << Bond("HEISENBERG", "J", {1, 0});
  bondlist << Bond("HEISENBERG", "J", {2, 3});
  bondlist << Bond("HEISENBERG", "J", {3, 2});
  bondlist << Bond("HEISENBERG", "J", {0, 2});
  bondlist << Bond("HEISENBERG", "J", {2, 0});
  bondlist << Bond("HEISENBERG", "J", {1, 3});
  bondlist << Bond("HEISENBERG", "J", {3, 1});
  return std::make_tuple(bondlist, couplings);
}

std::tuple<BondList, Couplings> square3x3(double t, double J) {
  BondList bondlist;
  Couplings couplings;
  couplings["T"] = t;
  couplings["J"] = J;
  bondlist << Bond("HUBBARDHOP", "T", {0, 1});
  bondlist << Bond("HUBBARDHOP", "T", {1, 2});
  bondlist << Bond("HUBBARDHOP", "T", {2, 0});
  bondlist << Bond("HUBBARDHOP", "T", {3, 4});
  bondlist << Bond("HUBBARDHOP", "T", {4, 5});
  bondlist << Bond("HUBBARDHOP", "T", {5, 3});
  bondlist << Bond("HUBBARDHOP", "T", {6, 7});
  bondlist << Bond("HUBBARDHOP", "T", {7, 8});
  bondlist << Bond("HUBBARDHOP", "T", {8, 6});
  bondlist << Bond("HUBBARDHOP", "T", {0, 3});
  bondlist << Bond("HUBBARDHOP", "T", {3, 6});
  bondlist << Bond("HUBBARDHOP", "T", {6, 0});
  bondlist << Bond("HUBBARDHOP", "T", {1, 4});
  bondlist << Bond("HUBBARDHOP", "T", {4, 7});
  bondlist << Bond("HUBBARDHOP", "T", {7, 1});
  bondlist << Bond("HUBBARDHOP", "T", {2, 5});
  bondlist << Bond("HUBBARDHOP", "T", {5, 8});
  bondlist << Bond("HUBBARDHOP", "T", {8, 2});
  bondlist << Bond("HEISENBERG", "J", {0, 1});
  bondlist << Bond("HEISENBERG", "J", {1, 2});
  bondlist << Bond("HEISENBERG", "J", {2, 0});
  bondlist << Bond("HEISENBERG", "J", {3, 4});
  bondlist << Bond("HEISENBERG", "J", {4, 5});
  bondlist << Bond("HEISENBERG", "J", {5, 3});
  bondlist << Bond("HEISENBERG", "J", {6, 7});
  bondlist << Bond("HEISENBERG", "J", {7, 8});
  bondlist << Bond("HEISENBERG", "J", {8, 6});
  bondlist << Bond("HEISENBERG", "J", {0, 3});
  bondlist << Bond("HEISENBERG", "J", {3, 6});
  bondlist << Bond("HEISENBERG", "J", {6, 0});
  bondlist << Bond("HEISENBERG", "J", {1, 4});
  bondlist << Bond("HEISENBERG", "J", {4, 7});
  bondlist << Bond("HEISENBERG", "J", {7, 1});
  bondlist << Bond("HEISENBERG", "J", {2, 5});
  bondlist << Bond("HEISENBERG", "J", {5, 8});
  bondlist << Bond("HEISENBERG", "J", {8, 2});
  return std::make_tuple(bondlist, couplings);
}

} // namespace hydra::testcases::electron
