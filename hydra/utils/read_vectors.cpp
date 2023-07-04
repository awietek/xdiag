#include "read_vectors.h"

#include <extern/fmt/format.h>
#include <hydra/common.h>

namespace hydra {

template <typename coeff_t>
arma::Mat<coeff_t> read_vectors(std::string type, std::string path_to_vecs,
                                int n) {
  using namespace arma;
  using namespace fmt;

  // get number of vecs if n=-1
  if (n == 0) {
    n = 0;
    while (file_exists(format("{}/{}_{}.arm", type, path_to_vecs, n))) {
      ++n;
    }
  } else if (n < 0) {
    Log.err(format("Invalid argument for number n of \"{}\" vectors", type));
  }

  if (n == 0) {
    return Mat<coeff_t>();
  } else {
    std::string filename = format("{}/{}_{}.arm", type, path_to_vecs, 0);
    if (!file_exists(filename)) {
      Log.err("Unable to read \"{}\" vector from file {}", type, filename);
    }
    Col<coeff_t> v;
    v.load(filename);
    int64_t N = v.size();

    Mat<coeff_t> Avecs(N, n);
    Avecs.col(0) = v;

    for (int i = 1; i < n; ++i) {
      std::string filename = format("{}/{}_{}.arm", type, path_to_vecs, n);
      if (!file_exists(filename)) {
        Log.err("Unable to read \"{}\" vector from file {}", type, filename);
      }
      v.load(filename);
      Avecs.col(n) = v;
    }

    return Avecs;
  }
}

template arma::Mat<double>
read_vectors<double>(std::string type, std::string path_to_vecs, int n);
template arma::Mat<complex>
read_vectors<complex>(std::string type, std::string path_to_vecs, int n);

} // namespace hydra
