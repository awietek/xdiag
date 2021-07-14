#include "stable_dot.h"

#include <mpi.h>
#include <cmath>
#include <cstring>

namespace hydra::mpi {

template <class coeff_t>
coeff_t StableDot(lila::Vector<coeff_t> const &v,
                  lila::Vector<coeff_t> const &w) {
  assert(v.size() == w.size());
  uint64 size = v.size();
  return detail::stable_dot_product(size, v.data(), w.data());
}

template float StableDot<float>(lila::Vector<float> const &v,
                                  lila::Vector<float> const &w);
template double StableDot<double>(lila::Vector<double> const &v,
                                  lila::Vector<double> const &w);
template scomplex StableDot<scomplex>(lila::Vector<scomplex> const &v,
                                  lila::Vector<scomplex> const &w);
template complex StableDot<complex>(lila::Vector<complex> const &v,
                                  lila::Vector<complex> const &w);




template <class coeff_t>
lila::real_t<coeff_t> StableNorm(lila::Vector<coeff_t> const &v) {
  return lila::real(std::sqrt(StableDot(v, v)));
}

template float StableNorm<float>(lila::Vector<float> const &v);
template double StableNorm<double>(lila::Vector<double> const &v);
template float StableNorm<scomplex>(lila::Vector<scomplex> const &v);
template double StableNorm<complex>(lila::Vector<complex> const &v);



} // namespace hydra::mpi

namespace hydra::mpi::detail {

// stable_dot_product does a dot product in parallel
// delivering always an identical result, no matter
// how the vectors are subdivided.

static inline void normalize_sums(int64 *isum) {
  for (int j = 0; j < 3; j++) {
    if (isum[j] < 0) {
      isum[j] += ((1ULL) << 60);
      isum[j + 1] -= ((1ULL) << 30);
    }
    isum[j + 1] += isum[j] >> 30;
    isum[j] &= 0x3fffffff;
  }
}

inline double stable_dot_product(const uint64 &n, const double *x,
                                 const double *y) {
  double absmax, tmp, fact, r_fact, sig;
  uint64 i;
  int64 itmp[4], isum[4], ival1, ival2;
  int iexp;

  absmax = 0.;
  for (i = 0; i < n; i++) {
    tmp = x[i] * y[i];
    if (fabs(tmp) > absmax)
      absmax = fabs(tmp);
  }

  tmp = absmax;
  MPI_Allreduce(&tmp, &absmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  frexp(absmax, &iexp); // get exponent of absmax

  fact = ldexp(1., 30 - iexp);   // same as 2**(30-iexp)
  r_fact = ldexp(1., iexp - 30); // 1./fact

  memset(isum, 0, sizeof(isum));

  for (i = 0; i < n; i++) {

    // Scale x[i]*y[i] into range -2**30 < tmp < 2**30
    // and store integer part in ival1
    tmp = (x[i] * y[i]) * fact;
    ival1 = tmp;

    // Scale fraction by 2**30 and store integer part in ival2

    ival2 = (tmp - ival1) * 1073741824.;

    isum[0] += ival2;
    isum[1] += ival1;

    // normalize sums from time to time (at least every 2**30 iterations)
    if ((n & 0xffff) == 0xffff)
      normalize_sums(isum);
  }

  normalize_sums(isum);
  memcpy(itmp, isum, sizeof(isum));
  MPI_Allreduce(&itmp, &isum, 4, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  normalize_sums(isum);

  sig = 1.;
  if (isum[3] < 0) {
    isum[3] = -isum[3];
    isum[2] = -isum[2];
    isum[1] = -isum[1];
    isum[0] = -isum[0];
    normalize_sums(isum);
    sig = -1.;
  }

  tmp = ldexp((double)isum[3], 60) + ldexp((double)isum[2], 30) +
        (double)isum[1] + ldexp((double)isum[0], -30);

  return sig * tmp * r_fact;
}

inline double stable_dot_product(const uint64 &n, const float *x,
                                 const float *y) {
  double absmax, tmp, fact, r_fact, sig;
  uint64 i;
  int64 itmp[4], isum[4], ival1, ival2;
  int iexp;

  absmax = 0.;
  for (i = 0; i < n; i++) {
    tmp = double(x[i]) * double(y[i]);

    if (fabs(tmp) > absmax)
      absmax = fabs(tmp);
  }

  tmp = absmax;
  MPI_Allreduce(&tmp, &absmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  frexp(absmax, &iexp); // get exponent of absmax

  fact = ldexp(1., 30 - iexp);   // same as 2**(30-iexp)
  r_fact = ldexp(1., iexp - 30); // 1./fact

  memset(isum, 0, sizeof(isum));

  for (i = 0; i < n; i++) {

    // Scale x[i]*y[i] into range -2**30 < tmp < 2**30
    // and store integer part in ival1

    tmp = double(x[i]) * double(y[i]) * fact;

    ival1 = tmp;

    // Scale fraction by 2**30 and store integer part in ival2

    ival2 = (tmp - ival1) * 1073741824.;

    isum[0] += ival2;
    isum[1] += ival1;

    // normalize sums from time to time (at least every 2**30 iterations)
    if ((n & 0xffff) == 0xffff)
      normalize_sums(isum);
  }

  normalize_sums(isum);
  memcpy(itmp, isum, sizeof(isum));
  MPI_Allreduce(&itmp, &isum, 4, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  normalize_sums(isum);

  sig = 1.;
  if (isum[3] < 0) {
    isum[3] = -isum[3];
    isum[2] = -isum[2];
    isum[1] = -isum[1];
    isum[0] = -isum[0];
    normalize_sums(isum);
    sig = -1.;
  }

  tmp = ldexp((double)isum[3], 60) + ldexp((double)isum[2], 30) +
        (double)isum[1] + ldexp((double)isum[0], -30);

  return sig * tmp * r_fact;
}

inline complex stable_dot_product(const uint64 &n, const complex *x,
                                  const complex *y) {
  double absmax, tmp, fact, r_fact, sig;
  uint64 i;
  int64 itmp[4], isum[4], ival1, ival2;
  int iexp;

  complex ctmp;

  absmax = 0.;
  for (i = 0; i < n; i++) {
    ctmp = std::conj(x[i]) * y[i];
    if (std::abs(std::real(ctmp)) > absmax)
      absmax = std::abs(std::real(ctmp));
    if (std::abs(std::imag(ctmp)) > absmax)
      absmax = std::abs(std::imag(ctmp));
  }

  tmp = absmax;
  MPI_Allreduce(&tmp, &absmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  frexp(absmax, &iexp); // get exponent of absmax

  fact = ldexp(1., 30 - iexp);   // same as 2**(30-iexp)
  r_fact = ldexp(1., iexp - 30); // 1./fact

  //
  // first do the real part
  //

  memset(isum, 0, sizeof(isum));

  for (i = 0; i < n; i++) {

    // Scale x[i]*y[i] into range -2**30 < tmp < 2**30
    // and store integer part in ival1

    tmp = std::real(std::conj(x[i]) * y[i]) * fact;
    ival1 = tmp;

    // Scale fraction by 2**30 and store integer part in ival2

    ival2 = (tmp - ival1) * 1073741824.;

    isum[0] += ival2;
    isum[1] += ival1;

    // normalize sums from time to time (at least every 2**30 iterations)
    if ((n & 0xffff) == 0xffff)
      normalize_sums(isum);
  }

  normalize_sums(isum);
  memcpy(itmp, isum, sizeof(isum));
  MPI_Allreduce(&itmp, &isum, 4, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  normalize_sums(isum);

  sig = 1.;
  if (isum[3] < 0) {
    isum[3] = -isum[3];
    isum[2] = -isum[2];
    isum[1] = -isum[1];
    isum[0] = -isum[0];
    normalize_sums(isum);
    sig = -1.;
  }

  tmp = ldexp((double)isum[3], 60) + ldexp((double)isum[2], 30) +
        (double)isum[1] + ldexp((double)isum[0], -30);

  double realpart = sig * tmp * r_fact;

  //
  // now do the imaginary part
  //

  memset(isum, 0, sizeof(isum));

  for (i = 0; i < n; i++) {

    // Scale x[i]*y[i] into range -2**30 < tmp < 2**30
    // and store integer part in ival1

    tmp = std::imag(std::conj(x[i]) * y[i]) * fact;
    ival1 = tmp;

    // Scale fraction by 2**30 and store integer part in ival2

    ival2 = (tmp - ival1) * 1073741824.;

    isum[0] += ival2;
    isum[1] += ival1;

    // normalize sums from time to time (at least every 2**30 iterations)
    if ((n & 0xffff) == 0xffff)
      normalize_sums(isum);
  }

  normalize_sums(isum);
  memcpy(itmp, isum, sizeof(isum));
  MPI_Allreduce(&itmp, &isum, 4, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  normalize_sums(isum);

  sig = 1.;
  if (isum[3] < 0) {
    isum[3] = -isum[3];
    isum[2] = -isum[2];
    isum[1] = -isum[1];
    isum[0] = -isum[0];
    normalize_sums(isum);
    sig = -1.;
  }

  tmp = ldexp((double)isum[3], 60) + ldexp((double)isum[2], 30) +
        (double)isum[1] + ldexp((double)isum[0], -30);

  double imagpart = sig * tmp * r_fact;
  return complex(realpart, imagpart);
}

inline scomplex stable_dot_product(const uint64 &n, const scomplex *x,
                                   const scomplex *y) {
  double absmax, tmp, fact, r_fact, sig;
  uint64 i;
  int64 itmp[4], isum[4], ival1, ival2;
  int iexp;

  complex ctmp;

  absmax = 0.;
  for (i = 0; i < n; i++) {
    ctmp = std::conj(x[i]) * y[i];
    if (std::abs(std::real(ctmp)) > absmax)
      absmax = std::abs(std::real(ctmp));
    if (std::abs(std::imag(ctmp)) > absmax)
      absmax = std::abs(std::imag(ctmp));
  }

  tmp = absmax;
  MPI_Allreduce(&tmp, &absmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  frexp(absmax, &iexp); // get exponent of absmax

  fact = ldexp(1., 30 - iexp);   // same as 2**(30-iexp)
  r_fact = ldexp(1., iexp - 30); // 1./fact

  //
  // first do the real part
  //

  memset(isum, 0, sizeof(isum));

  for (i = 0; i < n; i++) {

    // Scale x[i]*y[i] into range -2**30 < tmp < 2**30
    // and store integer part in ival1

    tmp = std::real(std::conj(x[i]) * y[i]) * fact;
    ival1 = tmp;

    // Scale fraction by 2**30 and store integer part in ival2

    ival2 = (tmp - ival1) * 1073741824.;

    isum[0] += ival2;
    isum[1] += ival1;

    // normalize sums from time to time (at least every 2**30 iterations)
    if ((n & 0xffff) == 0xffff)
      normalize_sums(isum);
  }

  normalize_sums(isum);
  memcpy(itmp, isum, sizeof(isum));
  MPI_Allreduce(&itmp, &isum, 4, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  normalize_sums(isum);

  sig = 1.;
  if (isum[3] < 0) {
    isum[3] = -isum[3];
    isum[2] = -isum[2];
    isum[1] = -isum[1];
    isum[0] = -isum[0];
    normalize_sums(isum);
    sig = -1.;
  }

  tmp = ldexp((double)isum[3], 60) + ldexp((double)isum[2], 30) +
        (double)isum[1] + ldexp((double)isum[0], -30);

  double realpart = sig * tmp * r_fact;

  //
  // now do the imaginary part
  //

  memset(isum, 0, sizeof(isum));

  for (i = 0; i < n; i++) {

    // Scale x[i]*y[i] into range -2**30 < tmp < 2**30
    // and store integer part in ival1

    tmp = std::imag(std::conj(x[i]) * y[i]) * fact;
    ival1 = tmp;

    // Scale fraction by 2**30 and store integer part in ival2

    ival2 = (tmp - ival1) * 1073741824.;

    isum[0] += ival2;
    isum[1] += ival1;

    // normalize sums from time to time (at least every 2**30 iterations)
    if ((n & 0xffff) == 0xffff)
      normalize_sums(isum);
  }

  normalize_sums(isum);
  memcpy(itmp, isum, sizeof(isum));
  MPI_Allreduce(&itmp, &isum, 4, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  normalize_sums(isum);

  sig = 1.;
  if (isum[3] < 0) {
    isum[3] = -isum[3];
    isum[2] = -isum[2];
    isum[1] = -isum[1];
    isum[0] = -isum[0];
    normalize_sums(isum);
    sig = -1.;
  }

  tmp = ldexp((double)isum[3], 60) + ldexp((double)isum[2], 30) +
        (double)isum[1] + ldexp((double)isum[0], -30);

  double imagpart = sig * tmp * r_fact;
  return scomplex(realpart, imagpart);
}

} // namespace hydra::mpi::detail
