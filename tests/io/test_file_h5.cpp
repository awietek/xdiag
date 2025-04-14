// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/common.hpp>
#include <xdiag/io/file_h5.hpp>

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

TEST_CASE("file_h5_append", "[io][hdf5]"){
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

TEST_CASE("fiile_h5_slicing", "[io][hdf5]"){
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

