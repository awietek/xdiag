// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "norm.hpp"

#include <xdiag/math/dot.hpp>
#include <xdiag/utils/error.hpp>

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/mpi/allreduce.hpp>
#endif

namespace xdiag::math {
template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v) try {
  return std::sqrt(xdiag::real(dot(block, v, v)));
}
XDIAG_CATCH

template double norm(Block const &, arma::Col<double> const &);
template double norm(Block const &, arma::Col<complex> const &);

template <typename coeff_t>
double norm1(Block const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, 1);
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return nrm;
}
XDIAG_CATCH

template double norm1(Block const &, arma::Col<double> const &);
template double norm1(Block const &, arma::Col<complex> const &);

template <typename coeff_t>
double norminf(Block const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, "inf");
#ifdef XDIAG_DISTRIBUTED
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  return nrm;
}
XDIAG_CATCH

template double norminf(Block const &, arma::Col<double> const &);
template double norminf(Block const &, arma::Col<complex> const &);

} // namespace xdiag::math
