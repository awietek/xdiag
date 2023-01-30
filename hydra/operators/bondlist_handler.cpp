#include "bondlist_handler.h"

#include <algorithm>
#include <hydra/utils/logger.h>

namespace hydra {

BondListHandler::BondListHandler(std::string key,
                                 std::map<std::string, complex> &couplings,
                                 std::map<std::string, arma::cx_mat> &matrices)
    : key_(key), couplings_(couplings), matrices_(matrices) {}

template <> double BondListHandler::as<double>() const {
  if (couplings_.find(key_) != couplings_.end()) {
    complex val = couplings_.at(key_);
    if (std::abs(val.imag()) < 1e-12) {
      return val.real();
    } else {
      Log.err("Error fetching real coupling from bond list: value is complex");
      return 0.;
    }
  } else {
    Log.err("Error fetching real coupling from bond list: key \"{}\" not found",
            key_);
    return 0.;
  }
}

template <> complex BondListHandler::as<complex>() const {
  if (couplings_.find(key_) != couplings_.end()) {
    return couplings_.at(key_);
  } else {
    Log.err("Error fetching coupling from bond list: key \"{}\" not found",
            key_);
    return 0.;
  }
}

template <> arma::mat BondListHandler::as<arma::mat>() const {
  if (matrices_.find(key_) != matrices_.end()) {
    arma::cx_mat matrix = matrices_.at(key_);
    if (arma::norm(arma::imag(matrix)) < 1e-12) {
      return arma::real(matrix);
    } else {
      Log.err("Error fetching real matrix from bond list: matrix is complex");
      return arma::mat();
    }
  } else {
    Log.err("Error fetching real matrix from bond list: key \"{}\" not found",
            key_);
    return arma::mat();
  }
}

template <> arma::cx_mat BondListHandler::as<arma::cx_mat>() const {
  if (matrices_.find(key_) != matrices_.end()) {
    return matrices_.at(key_);
  } else {
    Log.err("Error fetching matrix from bond list: key \"{}\" not found", key_);
    return arma::cx_mat();
  }
}

template <> void BondListHandler::operator=(int const &data) {
  couplings_[key_] = (double)data;
}
template <> void BondListHandler::operator=(double const &data) {
  couplings_[key_] = data;
}
template <> void BondListHandler::operator=(complex const &data) {
  couplings_[key_] = data;
}
template <> void BondListHandler::operator=(arma::mat const &data) {
  matrices_[key_] = arma::cx_mat(data, arma::zeros(data.n_rows, data.n_cols));
}
template <> void BondListHandler::operator=(arma::cx_mat const &data) {
  matrices_[key_] = data;
}

} // namespace hydra
