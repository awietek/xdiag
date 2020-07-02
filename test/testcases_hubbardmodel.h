#ifndef HYDRA_TEST_TESTCASES_HUBBARDMODEL
#define HYDRA_TEST_TESTCASES_HUBBARDMODEL

#include <hydra/all.h>
#include <random>

namespace hydra { namespace hubbardtestcases {
using namespace operators;
    
inline std::tuple<BondList, Couplings> heisenberg_triangle()
{
  BondList bondlist;
  bondlist << Bond("HEISENBERG", "J", {0, 1});
  bondlist << Bond("HEISENBERG", "J", {1, 2});
  bondlist << Bond("HEISENBERG", "J", {2, 0});

  Couplings couplings;
  couplings["J"] = 1.0;
  return std::make_tuple(bondlist, couplings);
}


inline std::tuple<BondList, Couplings> heisenberg_alltoall(int n_sites)
{
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  BondList bondlist;
  Couplings couplings;
  for (int s1=0; s1<n_sites; ++s1)
    for (int s2=s1+1; s2<n_sites; ++s2)
      {
	std::stringstream ss;
	ss << "J" << s1 << "_" << s2;
	std::string name = ss.str();
	double value = distribution(generator);
	bondlist << Bond("HEISENBERG", name, {s1, s2});
	couplings[name] = value;
      }
  return std::make_tuple(bondlist, couplings);
}

inline std::tuple<BondList, Couplings> heisenberg_kagome15()
{

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

inline std::tuple<BondList, Couplings> heisenberg_kagome39()
{
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



inline std::tuple<BondList, Couplings> freefermion_alltoall(int n_sites)
{
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  BondList bondlist;
  Couplings couplings;
  for (int s1=0; s1<n_sites; ++s1)
    for (int s2=s1+1; s2<n_sites; ++s2)
      {
	std::stringstream ss;
	ss << "T" << s1 << "_" << s2;
	std::string name = ss.str();
	double value = distribution(generator);
	bondlist << Bond("HUBBARDHOP", name, {s1, s2});
	couplings[name] = value;
      }
  return std::make_tuple(bondlist, couplings);
}


inline std::tuple<BondList, Couplings> freefermion_alltoall_complex(int n_sites)
{
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);

  BondList bondlist;
  Couplings couplings;
  for (int s1=0; s1<n_sites; ++s1)
    for (int s2=s1+1; s2<n_sites; ++s2)
      {
	std::stringstream ss;
	ss << "T" << s1 << "_" << s2;
	std::string name = ss.str();
	complex value = complex(distribution(generator),
				distribution(generator));
	bondlist << Bond("HUBBARDHOP", name, {s1, s2});
	couplings[name] = value;
      }
  return std::make_tuple(bondlist, couplings);
}

inline std::tuple<BondList, Couplings> tJchain(int n_sites, double t, double J)
{

  BondList bondlist;
  Couplings couplings;
  couplings["T"] = t;
  couplings["J"] = J;
  for (int s=0; s<n_sites; ++s)
    {
      bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % n_sites});
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    }
  return std::make_tuple(bondlist, couplings);
}


inline std::tuple<BondList, Couplings> randomAlltoAll4NoU()
{
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
  return std::make_tuple(bondlist, couplings);
}

inline std::tuple<BondList, Couplings> randomAlltoAll4()
{
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
  return std::make_tuple(bondlist, couplings);
}



inline std::tuple<BondList, Couplings> randomAlltoAll3()
{
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


inline std::tuple<BondList, Couplings> square2x2(double t, double J)
{
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

inline std::tuple<BondList, Couplings> square3x3(double t, double J)
{
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

    
}}
#endif
