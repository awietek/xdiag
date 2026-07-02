// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/config.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/io/file_h5.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>

#ifdef XDIAG_USE_HDF5
TEST_CASE("file_h5", "[io][hdf5]") {
  using namespace xdiag;
  using namespace arma;

  std::string filename = XDIAG_DIRECTORY "/misc/data/hdf5/write.h5";
  auto fl = FileH5(filename, "w!");
  fl["val"] = 12;
  fl["test/to"] = 22;
  fl["test/to2/group"] = 32;
  fl["test/to3/group2/asdf"] = 42;

  auto mat = arma::cx_mat(3, 5, arma::fill::randn);
  fl["a/b/c/mat"] = mat;
}

TEST_CASE("file_h5_append", "[io][hdf5]") {
  using namespace xdiag;
  using namespace arma;

  std::string filename = XDIAG_DIRECTORY "/misc/data/hdf5/write.h5";
  auto fid = FileH5(filename, "a");

  mat A(10, 5, arma::fill::zeros);
  A.col(1) = arma::vec(10, arma::fill::ones);
  arma::vec a(10, arma::fill::randn);

  fid["A"] = A;
  fid["A"].col(2) = a;
}

TEST_CASE("fiile_h5_slicing", "[io][hdf5]") {
  using namespace arma;
  using namespace xdiag;

  std::string filename = XDIAG_DIRECTORY "/misc/data/hdf5/write.h5";
  auto fid = FileH5(filename, "w!");

  arma::Cube<double> cubie(4, 5, 6, arma::fill::zeros);
  fid["cubie"] = cubie;

  auto sub_cube = arma::Mat<double>(4, 5, arma::fill::ones);
  fid["sub_cube"] = sub_cube;

  fid["cubie"].slice(3) = sub_cube;
}

#endif

#ifdef XDIAG_USE_HDF5
TEST_CASE("file_h5 scalar and vector writes", "[io][hdf5]") {
  using namespace xdiag;
  using namespace arma;

  std::string filename = XDIAG_DIRECTORY "/misc/data/hdf5/write_types.h5";
  auto fl = FileH5(filename, "w!");

  // --- scalar types (each triggers a distinct write_scalar instantiation) ---
  fl["s/i32"] = (int32_t)7;
  fl["s/i64"] = (int64_t)9;
  fl["s/u32"] = (uint32_t)11;
  fl["s/u64"] = (uint64_t)13;
  fl["s/dbl"] = 3.14;
  fl["s/cplx"] = complex(1.0, -2.0);

  // --- std::vector writes ---
  fl["v/ints"] = std::vector<int64_t>{1, 2, 3, 4};
  fl["v/dbls"] = std::vector<double>{1.5, 2.5, 3.5};
  fl["v/cplx"] = std::vector<complex>{{1.0, 2.0}, {3.0, 4.0}};

  // --- armadillo vector writes (real and complex) ---
  fl["a/vec"] = arma::vec{1.0, 2.0, 3.0};
  fl["a/cxvec"] = arma::cx_vec(arma::vec{1., 2.}, arma::vec{3., 4.});
  fl["a/cxmat"] = arma::cx_mat(2, 2, arma::fill::randn);
}

TEST_CASE("file_h5 overwrite existing dataset", "[io][hdf5]") {
  using namespace xdiag;
  using namespace arma;

  std::string filename = XDIAG_DIRECTORY "/misc/data/hdf5/write_overwrite.h5";
  auto fl = FileH5(filename, "w!");

  // Writing the same field twice hits the "dataset already exists" branch.
  fl["dup"] = 1.0;
  fl["dup"] = 2.0;

  arma::mat A(3, 3, arma::fill::ones);
  arma::mat B = 2.0 * A;
  fl["mat"] = A;
  fl["mat"] = B;
}
#endif
