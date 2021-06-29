#ifndef HYDRA_QNS_QN_TJ_H_
#define HYDRA_QNS_QN_TJ_H_

#include <string>
#include <sstream>

#include <hydra/common.h>

namespace hydra {

// quantum numbers "nup/ndn" for tj models
struct qn_tj {
  number_t n_up;
  number_t n_dn;
};

// Comparison operators
inline bool operator==(qn_tj const &q1, qn_tj const &q2) {
  return (q1.n_up == q2.n_up) && (q1.n_dn == q2.n_dn);
}
inline bool operator<(qn_tj const &q1, qn_tj const &q2) {
  return (q1.n_up < q2.n_up) || ((q1.n_up == q2.n_up) && (q1.n_dn < q2.n_dn));
}
inline bool operator!=(qn_tj const &q1, qn_tj const &q2) { return !(q1 == q2); }
inline bool operator<=(qn_tj const &q1, qn_tj const &q2) {
  return q1 < q2 || q1 == q2;
}
inline bool operator>(qn_tj const &q1, qn_tj const &q2) { return !(q1 <= q2); }
inline bool operator>=(qn_tj const &q1, qn_tj const &q2) { return !(q1 < q2); }

// Arithmetic operations
inline qn_tj operator+(qn_tj const &q1, qn_tj const &q2) {
  return qn_tj({q1.n_up + q2.n_up, q1.n_dn + q2.n_dn});
}
inline qn_tj operator-(qn_tj const &q1, qn_tj const &q2) {
  return qn_tj({q1.n_up - q2.n_up, q1.n_dn - q2.n_dn});
}
inline bool valid(qn_tj const &qn, number_t const &n_sites) {
  // (unsigned cast checks if positive)
  return ((number_t)qn.n_up <= n_sites) && ((number_t)qn.n_dn <= n_sites) &&
         (qn.n_up + qn.n_dn <= n_sites);
}

inline std::string String(qn_tj qn) {
  std::stringstream ss;
  ss << "QN TJ (n_up=" << qn.n_up << ", n_dn=" << qn.n_dn << ")";
  return ss.str();
}

} // namespace hydra

#endif
