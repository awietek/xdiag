#pragma once

#include "extern/armadillo/armadillo"

#include <string>

#include <hydra/algorithms/lanczos/tmatrix.h>

#include <hydra/blocks/blocks.h>
#include <hydra/blocks/electron/electron.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>

#include <hydra/operators/bondlist.h>
#include <hydra/states/product_state.h>
#include <hydra/states/random_state.h>
#include <hydra/states/state.h>
#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::utils {

void PrintPretty(const char *identifier, std::string str);

void PrintPretty(const char *identifier, int number);
void PrintPretty(const char *identifier, uint32_t number);
void PrintPretty(const char *identifier, uint64_t number);
void PrintPretty(const char *identifier, int64_t number);
void PrintPretty(const char *identifier, double number);
void PrintPretty(const char *identifier, complex number);

void PrintPretty(const char *identifier, Bond const &bond);
void PrintPretty(const char *identifier, BondList const &bondlist);
void PrintPretty(const char *identifier, Permutation const &perm);
void PrintPretty(const char *identifier, PermutationGroup const &group);
void PrintPretty(const char *identifier, Representation const &irrep);

void PrintPretty(const char *identifier, Block const &block);
void PrintPretty(const char *identifier, Spinhalf const &block);
void PrintPretty(const char *identifier, tJ const &block);
void PrintPretty(const char *identifier, Electron const &block);

void PrintPretty(const char *identifier, RandomState const &rstate);
void PrintPretty(const char *identifier, ProductState const &pstate);

void PrintPretty(const char *identifier, std::vector<float> const &vec);
void PrintPretty(const char *identifier, std::vector<double> const &vec);
void PrintPretty(const char *identifier, std::vector<scomplex> const &vec);
void PrintPretty(const char *identifier, std::vector<complex> const &vec);

void PrintPretty(const char *identifier, arma::Mat<float> const &mat);
void PrintPretty(const char *identifier, arma::Mat<double> const &mat);
void PrintPretty(const char *identifier,
                 arma::Mat<std::complex<float>> const &mat);
void PrintPretty(const char *identifier,
                 arma::Mat<std::complex<double>> const &mat);

void PrintPretty(const char *identifier, arma::Col<float> const &vec);
void PrintPretty(const char *identifier, arma::Col<double> const &vec);
void PrintPretty(const char *identifier,
                 arma::Col<std::complex<float>> const &vec);
void PrintPretty(const char *identifier,
                 arma::Col<std::complex<double>> const &vec);

void PrintPretty(const char *identifier, Tmatrix const &tmat);

template <typename coeff_t>
void PrintPretty(const char *identifier, State<coeff_t> const &state) {
  printf("%s:\n", identifier);
  if constexpr (is_real<coeff_t>()) {
    printf("  real state\n");
  } else {
    printf("  cplx state\n");
  }
  std::stringstream ss;
  ss.imbue(std::locale("en_US.UTF-8"));
  ss << state.size();
  printf("  dimension: %s\n", ss.str().c_str());
  printf("  norm     : %.9e\n", arma::norm(state.vector()));
  PrintPretty((std::string(identifier) + std::string(".block()")).c_str(),
              state.block());
}

} // namespace hydra::utils
