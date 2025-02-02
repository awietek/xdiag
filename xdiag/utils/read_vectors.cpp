#include "read_vectors.hpp"

#include <sstream>

#include <xdiag/common.hpp>

namespace xdiag {

static bool file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  return infile.good();
}

template <typename coeff_t>
arma::Mat<coeff_t> read_vectors(std::string type, std::string path_to_vecs,
                                int n) try {
  using namespace arma;
  using namespace fmt;

  std::stringstream ss;
  ss << type << "/" << path_to_vecs << "_" << n << ".arm";
  std::string filename = ss.str();

  // get number of vecs if n=-1
  if (n == 0) {
    n = 0;
    while (file_exists(filename)) {
      ++n;
    }
  } else if (n < 0) {
    XDIAG_THROW(
        format("Invalid argument for number n of \"{}\" vectors", type));
  }

  if (n == 0) {
    return Mat<coeff_t>();
  } else {
    if (!file_exists(filename)) {
      XDIAG_THROW(
          format("Unable to read \"{}\" vector from file {}", type, filename));
    }
    Col<coeff_t> v;
    v.load(filename);
    int64_t N = v.size();

    Mat<coeff_t> Avecs(N, n);
    Avecs.col(0) = v;

    for (int i = 1; i < n; ++i) {
      if (!file_exists(filename)) {
        XDIAG_THROW(format("Unable to read \"{}\" vector from file {}", type,
                           filename));
      }
      v.load(filename);
      Avecs.col(n) = v;
    }

    return Avecs;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template arma::Mat<double>
read_vectors<double>(std::string type, std::string path_to_vecs, int n);
template arma::Mat<complex>
read_vectors<complex>(std::string type, std::string path_to_vecs, int n);

} // namespace xdiag
