#include "arma_vector.hpp"

#include <xdiag/common.hpp>
#include <xdiag/io/toml/std_vector.hpp>
#include <xdiag/utils/type_string.hpp>

namespace xdiag::io {

template <typename T> arma::Col<T> arma_vector(toml::node const &node) try {
  return arma::Col<T>(std_vector<T>(node));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template arma::Col<double> arma_vector<double>(toml::node const &);
template arma::Col<complex> arma_vector<complex>(toml::node const &);
template arma::Col<arma::sword> arma_vector<arma::sword>(toml::node const &);
template arma::Col<arma::uword> arma_vector<arma::uword>(toml::node const &);

template <typename T> toml::array toml_array(arma::Col<T> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(x);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> toml::array toml_array(arma::Col<arma::uword> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back((int64_t)x);
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> toml::array toml_array(arma::Col<complex> const &value) try {
  toml::array arr;
  for (auto &&x : value) {
    arr.push_back(toml::array{x.real(), x.imag()});
  }
  return arr;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template toml::array toml_array(arma::Col<double> const &);
template toml::array toml_array(arma::Col<arma::sword> const &);

} // namespace xdiag::io
