#pragma once

#include <extern/armadillo/armadillo>

#include <string>
#include <map>
#include <hydra/common.h>

namespace hydra {
class BondListHandler {
public:
  BondListHandler(std::string key, std::map<std::string, complex> &couplings,
                  std::map<std::string, arma::cx_mat> &matrices);
  BondListHandler(BondListHandler const &) = delete;
  BondListHandler &operator=(BondListHandler const &) = delete;

  template <class data_t> data_t as() const;
  template <class data_t> void operator=(data_t const &data);

private:
  std::string key_;
  std::map<std::string, complex> &couplings_;
  std::map<std::string, arma::cx_mat> &matrices_;
};

} // namespace hydra
