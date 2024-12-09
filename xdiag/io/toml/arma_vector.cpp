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

} // namespace xdiag::io::toml
