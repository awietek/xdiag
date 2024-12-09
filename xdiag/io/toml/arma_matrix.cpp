#include "arma_vector.hpp"

#include <xdiag/common.hpp>
#include <xdiag/io/toml/std_vector.hpp>
#include <xdiag/io/toml/value.hpp>
#include <xdiag/utils/type_string.hpp>

namespace xdiag::io {

template <typename T> arma::Mat<T> arma_matrix(toml::node const &node) try {
  auto array = node.as_array();
  if (array) {
    toml::array rows = *array;
    std::size_t m = rows.size();
    std::size_t n = 0;
    if (m > 0) {
      std::vector<T> vector;
      for (std::size_t i = 0; i < m; ++i) {
        auto row_entries_opt = rows[i].as_array();
        if (row_entries_opt) {
          toml::array row_entries = *row_entries_opt;
          if (i == 0) {
            n = row_entries.size();
            vector.resize(m * n);
          }
          if (row_entries.size() != n) {
            XDIAG_THROW(
                fmt::format("Error reading TOML array to \"{}\": not all rows "
                            "have same length",
                            utils::type_string<arma::Mat<T>>()));
          }

          for (std::size_t j = 0; j < n; ++j) {
            vector[i + m * j] = value<T>(row_entries[j]);
          }

        } else {
          XDIAG_THROW(fmt::format(
              "Error reading TOML array to \"{}\": Cannot parse row {}!",
              utils::type_string<arma::Mat<T>>(), i));
        }
      } // for (int i = 0; i < m; ++i)
      return arma::Mat<T>(vector.data(), m, n);
    } else {
      return arma::Mat<T>();
    }
  } else {
    XDIAG_THROW(fmt::format(
        "TOML node cannot be converted to \"{}\". Node is not an array.",
        utils::type_string<arma::Mat<T>>()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template arma::Mat<double> arma_matrix<double>(toml::node const &);
template arma::Mat<complex> arma_matrix<complex>(toml::node const &);
template arma::Mat<arma::sword> arma_matrix<arma::sword>(toml::node const &);
template arma::Mat<arma::uword> arma_matrix<arma::uword>(toml::node const &);

} // namespace xdiag::io
