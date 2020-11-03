#ifndef HYDRA_QNS_QN_SPINHALF_H_
#define HYDRA_QNS_QN_SPINHALF_H_

#include <string>
#include <sstream>

#include <hydra/common.h>

namespace hydra {

// quantum numbers "nup/ndn" for spinhalf models
struct qn_spinhalf {
  number_t n_up;
};

// Comparison operators
inline bool operator==(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return (q1.n_up == q2.n_up);
}
inline bool operator<(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return (q1.n_up < q2.n_up);
}
inline bool operator!=(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return !(q1 == q2);
}
inline bool operator<=(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return q1 < q2 || q1 == q2;
}
inline bool operator>(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return !(q1 <= q2);
}
inline bool operator>=(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return !(q1 < q2);
}

// Arithmetic operations
inline qn_spinhalf operator+(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return qn_spinhalf({q1.n_up + q2.n_up});
}
inline qn_spinhalf operator-(qn_spinhalf const &q1, qn_spinhalf const &q2) {
  return qn_spinhalf({q1.n_up - q2.n_up});
}
inline bool valid(qn_spinhalf const &qn, number_t const &n_sites) {
  // (unsigned cast checks if positive)
  return ((unsigned)qn.n_up <= n_sites);
}
  
inline std::string String(qn_spinhalf qn) {
  std::stringstream ss;
  ss << "QN SpinHalf (n_up=" << qn.n_up << ")";
  return ss.str();
}
  
} // namespace hydra

#endif
