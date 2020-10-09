#ifndef HYDRA_QNS_QN_ELECTRON_H_
#define HYDRA_QNS_QN_ELECTRON_H_

#include <hydra/common.h>

namespace hydra {

// quantum numbers "nup/ndn" for electron models
struct qn_electron {
  number_t n_up;
  number_t n_dn;
};

// Comparison operators
inline bool operator==(qn_electron const &q1, qn_electron const &q2) {
  return (q1.n_up == q2.n_up) && (q1.n_dn == q2.n_dn);
}
inline bool operator<(qn_electron const &q1, qn_electron const &q2) {
  return (q1.n_up < q2.n_up) || ((q1.n_up == q2.n_up) && (q1.n_dn < q2.n_dn));
}
inline bool operator!=(qn_electron const &q1, qn_electron const &q2) {
  return !(q1 == q2);
}
inline bool operator<=(qn_electron const &q1, qn_electron const &q2) {
  return q1 < q2 || q1 == q2;
}
inline bool operator>(qn_electron const &q1, qn_electron const &q2) {
  return !(q1 <= q2);
}
inline bool operator>=(qn_electron const &q1, qn_electron const &q2) {
  return !(q1 < q2);
}

// Arithmetic operations
inline qn_electron operator+(qn_electron const &q1, qn_electron const &q2) {
  return qn_electron({q1.n_up + q2.n_up, q1.n_dn + q2.n_dn});
}
inline qn_electron operator-(qn_electron const &q1, qn_electron const &q2) {
  return qn_electron({q1.n_up - q2.n_up, q1.n_dn - q2.n_dn});
}
inline bool valid(qn_electron const &qn, number_t const &n_sites) {
  // (unsigned cast checks if positive)
  return ((unsigned)qn.n_up <= n_sites) && ((unsigned)qn.n_dn <= n_sites);
}

}

#endif
