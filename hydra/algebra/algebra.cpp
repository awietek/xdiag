#include "algebra.h"

#include <hydra/algebra/apply.h>

namespace hydra {

double norm(State const &v) try {
  if (v.isreal()) {
    return arma::norm(v.matrix(false));
  } else {
    return arma::norm(v.matrixC(false));
  }
} catch (...) {
  rethrow(__func__, "Unable compute norm of State");
  return 0.;
}

double dot(State const &v, State const &w) try {
  if ((v.isreal()) && (w.isreal())) {
    return arma::dot(v.matrix(false), w.matrix(false));
  } else {
    throw std::logic_error("Unable to compute real dot product of a complex "
                           "state, consider using dotC instead.");
  }
} catch (...) {
  rethrow(__func__, "Unable compute dot product of two States");
  return 0.;
}

complex dotC(State const &v, State const &w) try {
  if ((v.isreal()) && (w.isreal())) {
    return arma::dot(v.matrix(false), w.matrix(false));
  } else if ((v.isreal()) && (w.iscomplex())) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      throw(std::runtime_error("Unable to create intermediate complex State"));
    }
    return arma::cdot(v2.matrixC(false), w.matrixC(false));
  } else if ((w.isreal()) && (v.iscomplex())) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      throw(std::runtime_error("Unable to create intermediate complex State"));
    }
    return arma::cdot(v.matrixC(false), w.matrixC(false));
  } else {
    return arma::cdot(v.matrixC(false), w.matrixC(false));
  }
} catch (...) {
  rethrow(__func__, "Unable compute dot product of two States");
  return 0.;
}

double inner(BondList const &bonds, State const &v) try {
  auto w = v;
  apply(bonds, v, w);
  return dot(w, v);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"inner\"");
  return 0.;
}

complex innerC(BondList const &bonds, State const &v) try {
  auto w = v;
  apply(bonds, v, w);
  return dotC(w, v);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"innerC\"");
  return 0.;
}

double inner(Bond const &bond, State const &v) try {
  return inner(BondList({bond}), v);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"inner\"");
  return 0.;
}

complex innerC(Bond const &bond, State const &v) try {
  return innerC(BondList({bond}), v);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"innerC\"");
  return 0.;
}

double inner(State const &v, BondList const &bonds, State const &w) try {
  auto bw = zeros_like(w);
  apply(bonds, w, bw);
  return dot(v, bw);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"inner\"");
  return 0.;
}

complex innerC(State const &v, BondList const &bonds, State const &w) try {
  auto bw = zeros_like(w);
  apply(bonds, w, bw);
  return dotC(v, bw);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"innerC\"");
  return 0.;
}

double inner(State const &v, Bond const &bond, State const &w) try {
  return inner(v, BondList({bond}), w);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"inner\"");
  return 0.;
}

complex innerC(State const &v, Bond const &bond, State const &w) try {
  return innerC(v, BondList({bond}), w);
} catch (...) {
  rethrow(__func__, "Error computing expectation value using \"innerC\"");
  return 0.;
}

State &operator*=(State &X, complex alpha) {
  if (X.isreal()) {
    X.make_complex();
    X.matrixC(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
}
State &operator*=(State &X, double alpha) {
  if (X.isreal()) {
    X.matrix(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
}
State &operator/=(State &X, complex alpha) {
  if (X.isreal()) {
    X.make_complex();
    X.matrixC(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
}
State &operator/=(State &X, double alpha) {
  if (X.isreal()) {
    X.matrix(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
}

// template <class coeff_t>
// coeff_t dot(State<coeff_t> const &v, State<coeff_t> const &w) {
//   if (v.block() != w.block()) {
//     Log.err("Error: trying to perform Dot product of distinct blocks");
//   }
//   if constexpr (mpi::is_mpi_block<Block>) {
//     return DotMPI(v.vector(), w.vector());
//   } else {
//     return arma::cdot(v.vector(), w.vector());
//   }
// }
// template double dot(State<double> const &, State<double> const &);
// template complex dot(State<complex> const &, State<complex> const &);

// template <class coeff_t> real_t<coeff_t> norm(State<coeff_t> const &v) {
//   return std::abs(std::sqrt(dot(v, v)));
// }
// template double norm(State<double> const &);
// template double norm(State<complex> const &);

// template <class coeff_t>
// coeff_t inner(BondList const &bonds, State<coeff_t> const &v) {
//   auto Hv = State<coeff_t>(v.block());
//   apply(bonds, v, Hv);
//   return dot(v, Hv);
// }
// template double inner(BondList const &, State<double> const &);
// template complex inner(BondList const &, State<complex> const &);

// template <class coeff_t>
// coeff_t inner(Bond const &bond, State<coeff_t> const &v) {
//   BondList bonds;
//   bonds << bond;
//   return inner(bonds, v);
// }
// template double inner(Bond const &, State<double> const &);
// template complex inner(Bond const &, State<complex> const &);

// template <class coeff_t>
// coeff_t inner(State<coeff_t> const &w, BondList const &bonds,
//               State<coeff_t> const &v) {
//   auto Hv = State<coeff_t>(w.block());
//   apply(bonds, v, Hv);
//   return dot(w, Hv);
// }
// template double inner(State<double> const &, BondList const &,
//                       State<double> const &);
// template complex inner(State<complex> const &, BondList const &,
//                        State<complex> const &);

// template <class coeff_t>
// coeff_t inner(State<coeff_t> const &w, Bond const &bond,
//               State<coeff_t> const &v) {
//   BondList bonds;
//   bonds << bond;
//   return inner(w, bonds, v);
// }
// template double inner(State<double> const &, Bond const &,
//                       State<double> const &);
// template complex inner(State<complex> const &, Bond const &,
//                        State<complex> const &);

// template <class coeff_t>
// void apply(BondList const &bonds, Block const &block_in,
//            arma::Col<coeff_t> const &vec_in, Block const &block_out,
//            arma::Col<coeff_t> &vec_out) {

//   std::visit(
//       variant::overloaded{
//           [&bonds, &vec_in, &vec_out](Spinhalf const &blk_in,
//                                       Spinhalf const &blk_out) {
//             apply(bonds, blk_in, vec_in, blk_out, vec_out);
//           },
//           [&bonds, &vec_in, &vec_out](tJ const &blk_in, tJ const &blk_out) {
//             apply(bonds, blk_in, vec_in, blk_out, vec_out);
//           },
//           [&bonds, &vec_in, &vec_out](Electron const &blk_in,
//                                       Electron const &blk_out) {
//             apply(bonds, blk_in, vec_in, blk_out, vec_out);
//           },
//           [](auto const &blk_in, auto const &blk_out) {
//             Log.err("Error in apply: Invalid blocks or combination of
//             blocks"); (void)blk_in; (void)blk_out;
//           },
//       },
//       block_in.variant(), block_out.variant());
// }
// template void apply(BondList const &, Block const &, arma::Col<double> const
// &,
//                     Block const &, arma::Col<double> &);
// template void apply(BondList const &, Block const &, arma::Col<complex> const
// &,
//                     Block const &, arma::Col<complex> &);

// template <class coeff_t>
// void apply(BondList const &bonds, State<coeff_t> const &state_in,
//            State<coeff_t> &state_out) {
//   apply(bonds, state_in.block(), state_in.vector(), state_out.block(),
//         state_out.vector());
// }
// template void apply(BondList const &, State<double> const &, State<double>
// &); template void apply(BondList const &, State<complex> const &,
// State<complex> &);

// template <class coeff_t>
// void apply(Bond const &bond, State<coeff_t> const &state_in,
//            State<coeff_t> &state_out) {
//   BondList bonds;
//   bonds << bond;
//   apply(bonds, state_in, state_out);
// }
// template void apply(Bond const &, State<double> const &, State<double> &);
// template void apply(Bond const &, State<complex> const &, State<complex> &);

// // template <class coeff_t>
// // void apply(Bond const &bond, complex coeff, State<coeff_t> const
// &state_in,
// //            State<coeff_t> &state_out) {
// //   BondList bonds;
// //   bonds << bond;
// //   apply(bonds, state_in, state_out);
// // }

// // template <class coeff_t>
// // State<coeff_t> apply(Bond const &bond, State<coeff_t> const &state_in) {
// //   BondList bonds;
// //   bonds << bond;
// //   return apply(bonds, state_in);
// // }

// State<complex> &operator/=(State<complex> &X, complex alpha) {
//   X.vector() /= alpha;
//   return X;
// }

// State<complex> &operator/=(State<complex> &X, double alpha) {
//   X.vector() /= alpha;
//   return X;
// }

// State<double> &operator/=(State<double> &X, double alpha) {
//   X.vector() /= alpha;
//   return X;
// }

// State<complex> &operator*=(State<complex> &X, complex alpha) {
//   X.vector() *= alpha;
//   return X;
// }

// State<complex> &operator*=(State<complex> &X, double alpha) {
//   X.vector() *= alpha;
//   return X;
// }

// State<double> &operator*=(State<double> &X, double alpha) {
//   X.vector() *= alpha;
//   return X;
// }

} // namespace hydra
