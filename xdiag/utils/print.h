#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <string>

#include <xdiag/algorithms/lanczos/tmatrix.h>

#include <xdiag/blocks/blocks.h>
#include <xdiag/blocks/electron/electron.h>
#include <xdiag/blocks/spinhalf/spinhalf.h>
#include <xdiag/blocks/tj/tj.h>

#include <xdiag/operators/bondlist.h>
#include <xdiag/states/product_state.h>
#include <xdiag/states/random_state.h>
#include <xdiag/states/state.h>

#include <xdiag/symmetries/continuous_group.h>
#include <xdiag/symmetries/permutation.h>
#include <xdiag/symmetries/permutation_group.h>
#include <xdiag/symmetries/qn.h>
#include <xdiag/symmetries/representation.h>

namespace xdiag::utils {

void print_pretty(const char *identifier, std::string str);

void print_pretty(const char *identifier, int number);
void print_pretty(const char *identifier, uint32_t number);
void print_pretty(const char *identifier, uint64_t number);
void print_pretty(const char *identifier, int64_t number);
void print_pretty(const char *identifier, double number);
void print_pretty(const char *identifier, complex number);

void print_pretty(const char *identifier, Bond const &bond);
void print_pretty(const char *identifier, BondList const &bondlist);
void print_pretty(const char *identifier, Permutation const &perm);
void print_pretty(const char *identifier, PermutationGroup const &group);
void print_pretty(const char *identifier, Representation const &irrep);
void print_pretty(const char *identifier, U1 const &qn);
void print_pretty(const char *identifier, QNum const &qn);
void print_pretty(const char *identifier, QN const &qn);

void print_pretty(const char *identifier, block_variant_t const &block);
void print_pretty(const char *identifier, Spinhalf const &block);
void print_pretty(const char *identifier, tJ const &block);
void print_pretty(const char *identifier, Electron const &block);

#ifdef XDIAG_USE_MPI
void print_pretty(const char *identifier, tJDistributed const &block);
#endif

void print_pretty(const char *identifier, RandomState const &rstate);
void print_pretty(const char *identifier, ProductState const &pstate);

void print_pretty(const char *identifier, std::vector<int> const &vec);
void print_pretty(const char *identifier, std::vector<float> const &vec);
void print_pretty(const char *identifier, std::vector<double> const &vec);
void print_pretty(const char *identifier, std::vector<complex> const &vec);
void print_pretty(const char *identifier, std::vector<std::string> const &vec);

void print_pretty(const char *identifier, arma::uvec const &vec);
void print_pretty(const char *identifier, arma::ivec const &vec);
void print_pretty(const char *identifier, arma::umat const &mat);
void print_pretty(const char *identifier, arma::imat const &mat);

void print_pretty(const char *identifier, arma::Mat<float> const &mat);
void print_pretty(const char *identifier, arma::Mat<double> const &mat);
void print_pretty(const char *identifier,
                  arma::Mat<std::complex<float>> const &mat);
void print_pretty(const char *identifier,
                  arma::Mat<std::complex<double>> const &mat);

void print_pretty(const char *identifier, arma::Col<float> const &vec);
void print_pretty(const char *identifier, arma::Col<double> const &vec);
void print_pretty(const char *identifier,
                  arma::Col<std::complex<float>> const &vec);
void print_pretty(const char *identifier,
                  arma::Col<std::complex<double>> const &vec);

void print_pretty(const char *identifier, Tmatrix const &tmat);
void print_pretty(const char *identifier, State const &state);

} // namespace xdiag::utils
