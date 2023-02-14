#include "../catch.hpp"

#include <hydra/all.h>

template <typename T> void test_write_read(T val) {
  using namespace hydra;

  // Test writing
  std::string filename = HYDRA_DIRECTORY"/misc/data/toml/write.toml";
  auto fl = FileToml(filename, 'w');
  std::string key = "val";

  fl[key] = val;
  fl.close();

  fl = FileToml(filename, 'r');
  T val_r = fl[key].as<T>();

  // HydraPrint(val);
  // HydraPrint(val_r);
  REQUIRE(val == val_r);
}

TEST_CASE("file_toml", "[io]") {
  using namespace hydra;
  using namespace arma;

  // Just try to parse everything in the example toml file
  std::string filename = HYDRA_DIRECTORY"/misc/data/toml/read.toml";
  auto fl = FileToml(filename, 'r');

  // Parse String
  std::string n = "title";
  // Log("A");
  // Log("{}", fl[n].as<std::string>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<std::string>() == "TOML Beispiel");

  // Parse integers
  n = "datenbank.verbindungen_max";
  // Log("B");
  // Log("{}", fl[n].as<int32_t>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<int8_t>() == 42);
  REQUIRE(fl[n].as<int16_t>() == 42);
  REQUIRE(fl[n].as<int32_t>() == 42);
  REQUIRE(fl[n].as<int64_t>() == 42);
  REQUIRE(fl[n].as<uint8_t>() == 42);
  REQUIRE(fl[n].as<uint16_t>() == 42);
  REQUIRE(fl[n].as<uint32_t>() == 42);
  REQUIRE(fl[n].as<uint64_t>() == 42);

  // Parse bool
  n = "datenbank.aktiviert";
  // Log("{}", fl[n].as<bool>());
  // std::cout << fl[n].as<bool>() << "\n";
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<bool>() == true);

  n = "floating.pi";
  // Log("{}", fl[n].as<double>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<double>() == 3.1415);

  n = "floating.cplx";
  // Log("{}", fl[n].as<complex>());
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<complex>() == 2.3122 + .1237i);

  n = "vectors.ints";
  // for (auto i : fl[n].as<std::vector<int>>()) {
  //   Log("i: {}", i);
  // }
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<std::vector<int>>() == std::vector<int>{1, 3, 7, 4, 5});

  n = "vectors.floats";
  // for (auto i : fl[n].as<std::vector<double>>()) {
  //   Log("d: {}", i);
  // }
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<std::vector<double>>() ==
          std::vector<double>{1.234, 7.654, 3.1415});

  n = "vectors.cplxs";
  // for (auto i : fl[n].as<std::vector<complex>>()) {
  //   Log("d: {}", i);
  // }
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<std::vector<complex>>() ==
          std::vector<complex>{2.3122 + 0.1237i, -12.3122 + 20.1237i,
                               32.3122 - 40.1237i});

  n = "clients.hosts";
  // for (auto i : fl[n].as<std::vector<std::string>>()) {
  //   Log("s: {}", i);
  // }
  REQUIRE(fl.defined(n));
  REQUIRE(fl[n].as<std::vector<std::string>>() ==
          std::vector<std::string>{"alpha", "omega"});

  ivec iv{8001, 8001, 8002};
  uvec uv{8001, 8001, 8002};
  REQUIRE(all(fl["datenbank.ports"].as<ivec>() == iv));
  REQUIRE(all(fl["datenbank.ports"].as<uvec>() == uv));

  auto sx = fl["pauli.sx"].as<mat>();
  REQUIRE(close(sx, mat{{0.0, 0.5}, {0.5, 0.0}}));
  // HydraPrint(sx);

  auto sy = fl["pauli.sy"].as<cx_mat>();
  cx_mat syy{{0. + 0.i, 0. - 0.5i}, {0. + 0.5i, 0.0}};
  // HydraPrint(sy);
  // HydraPrint(syy);

  REQUIRE(close(sy, syy));

  auto sz = fl["pauli.sz"].as<mat>();
  REQUIRE(close(sz, mat{{0.5, 0.}, {0., -0.5}}));
  // HydraPrint(sz);

  auto szc = fl["pauli.sz"].as<cx_mat>();
  REQUIRE(close(szc, cx_mat{{0.5 + 0i, 0. + 0i}, {0. + 0i, -0.5 + 0i}}));
  // HydraPrint(szc);

  auto other = fl["pauli.other"].as<mat>();
  REQUIRE(close(other, mat{{1.0, 2.0, 3.0}, {5.0, 6.0, 7.0}}));
  // HydraPrint(other);

  int a = 42;
  double b = 1.234;
  complex c = 1.234 + 4.321i;
  std::vector<int> d{1, 5, 4, 3, 9, 123};
  std::vector<double> e{1.23, 5.43, 4.54, 3.65, 9.76, 12.3};
  std::vector<complex> f{1.23 + 5.43i, 4.54 + 3.65i, 9.76 + 12.3i};
  vec g{1.23, 5.43, 4.54, 3.65, 9.76, 12.3};
  cx_vec h{1.23 + 5.43i, 4.54 + 3.65i, 9.76 + 12.3i};
  std::vector<std::string> i{"Stoan", "Mario", "Chrissi", "Alex"};
  mat j = mat{{0.0, 0.5}, {0.5, 0.0}};
  cx_mat k = cx_mat{{0. + 0.i, 0. - 0.5i}, {0. + 0i, 0.0 - 0.5i}};
  imat l = imat{{1, 2, 3}, {5, 6, 7}};
  umat m = umat{{1, 2, 3}, {5, 6, 7}};

  test_write_read(a);
  test_write_read(b);
  test_write_read(c);
  test_write_read(d);
  test_write_read(e);
  test_write_read(f);
  test_write_read(i);

  filename = HYDRA_DIRECTORY"/misc/data/toml/write.toml";
  fl = FileToml(filename, 'w');
  fl["a"] = a;
  fl["b"] = b;
  fl["c"] = c;
  fl["d"] = d;
  fl["e"] = e;
  fl["f"] = f;
  fl["g"] = g;
  fl["h"] = h;
  fl["i"] = i;
  fl["j"] = j;
  fl["k"] = k;
  fl["l"] = l;
  fl["m"] = m;
  fl.close();

  fl = FileToml(filename, 'r');
  auto ar = fl["a"].as<int>();
  auto br = fl["b"].as<double>();
  auto cr = fl["c"].as<complex>();
  auto dr = fl["d"].as<std::vector<int>>();
  auto er = fl["e"].as<std::vector<double>>();
  auto fr = fl["f"].as<std::vector<complex>>();
  auto gr = fl["g"].as<vec>();
  auto hr = fl["h"].as<cx_vec>();
  auto ir = fl["i"].as<std::vector<std::string>>();
  auto jr = fl["j"].as<mat>();
  auto kr = fl["k"].as<cx_mat>();
  auto lr = fl["l"].as<imat>();
  auto mr = fl["m"].as<umat>();

  // HydraPrint(lr);
  REQUIRE(a == ar);
  REQUIRE(b == br);
  REQUIRE(c == cr);
  REQUIRE(d == dr);
  REQUIRE(e == er);
  REQUIRE(f == fr);
  REQUIRE(close(g, gr));
  REQUIRE(close(h, hr));
  REQUIRE(i == ir);
  REQUIRE(close(j, jr));
  REQUIRE(close(k, kr));

  auto p = Permutation({5, 3, 4, 2, 1, 0, 7, 6});
  fl["perm"] = p;
  test_write_read(p);

  std::string lfile = HYDRA_DIRECTORY"/misc/data/triangular.j1j2jch/"
                      "triangular.12.j1j2jch.sublattices.fsl.lat";
  auto group = PermutationGroup(read_permutations(lfile));
  test_write_read(group);

  auto irrep = read_representation(lfile, "X.C1.A");
  test_write_read(irrep);

  auto bond = Bond("HB", "H", 1);
  test_write_read(bond);

  bond = Bond("HB", "J1", {1, 2});
  test_write_read(bond);

  bond = Bond("HB", 1.0, 1);
  test_write_read(bond);

  bond = Bond("HB", 1.2, {1, 2});
  test_write_read(bond);

  bond = Bond("HB", 2.1 + 1.2i, 1);
  test_write_read(bond);

  bond = Bond("HB", 2.1 + 1.2i, {1, 2});
  test_write_read(bond);

  auto matr = arma::mat(3, 3, arma::fill::randu);

  bond = Bond(matr, "H", 1);
  test_write_read(bond);

  bond = Bond(matr, "J1", {1, 2});
  test_write_read(bond);

  bond = Bond(matr, 1.0, 1);
  test_write_read(bond);

  bond = Bond(matr, 1.2, {1, 2});
  test_write_read(bond);

  bond = Bond(matr, 2.1 + 1.2i, 1);
  test_write_read(bond);

  bond = Bond(matr, 2.1 + 1.2i, {1, 2});
  test_write_read(bond);

  auto matc = arma::cx_mat(3, 3, arma::fill::randu);

  bond = Bond(matc, "H", 1);
  test_write_read(bond);

  bond = Bond(matc, "J1", {1, 2});
  test_write_read(bond);

  bond = Bond(matc, 1.0, 1);
  test_write_read(bond);

  bond = Bond(matc, 1.2, {1, 2});
  test_write_read(bond);

  bond = Bond(matc, 2.1 + 1.2i, 1);
  test_write_read(bond);

  bond = Bond(matc, 2.1 + 1.2i, {1, 2});
  test_write_read(bond);

  auto bonds = read_bondlist(lfile);
  test_write_read(bonds);

  bonds["J1"] = 1.0;
  bonds["J2"] = 0.2 - 0.1i;
  test_write_read(bonds);

  bonds["M1"] = arma::mat(3, 4, arma::fill::randu);
  bonds["M2"] = arma::cx_mat(3, 4, arma::fill::randu);
  test_write_read(bonds);
  
}
