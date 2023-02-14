#include "../catch.hpp"

#include <hydra/all.h>

#ifdef HYDRA_USE_HDF5
TEST_CASE("file_h5", "[io]") {
  using namespace hydra;
  using namespace arma;

  std::string filename = HYDRA_DIRECTORY "/misc/data/hdf5/write.h5";
  auto fl = FileH5(filename, "w!");
  fl["val"] = 12;
  fl["test/to"] = 22;
  fl["test/to2/group"] = 32;
  fl["test/to3/group2/asdf"] = 42;

  auto mat = arma::cx_mat(3, 5, arma::fill::randn);
  HydraPrint(mat);

  fl["a/b/c/mat"] = mat;
}
#endif
