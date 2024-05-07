#include "algebra.hpp"

#include <xdiag/algebra/apply.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/allreduce.hpp>
#include <xdiag/parallel/mpi/cdot_distributed.hpp>
#endif

namespace xdiag {

double norm(State const &v) try {
  if (v.n_cols() > 1) {
    XDIAG_THROW("Cannot compute norm of state with more than one column");
    return 0;
  } else {
    if (v.isreal()) {
      return norm(v.block(), v.vector(0, false));
    } else {
      return norm(v.block(), v.vectorC(0, false));
    }
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double dot(State const &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
    return 0;
  }

  if ((v.n_cols() > 1) || (w.n_cols() > 1)) {
    XDIAG_THROW(
        "Cannot compute dot product of state with more than one column");
    return 0;
  }

  if ((v.isreal()) && (w.isreal())) {
    return dot(v.vector(0, false), w.vector(0, false));
  } else {
    XDIAG_THROW("Unable to compute real dot product of a complex "
                "state, consider using dotC instead.");
    return 0;
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

complex dotC(State const &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot form dot product for states on different blocks");
  }

  if ((v.n_cols() > 1) || (w.n_cols() > 1)) {
    XDIAG_THROW(
        "Cannot compute dotC product of state with more than one column");
  }
  if ((v.isreal()) && (w.isreal())) {
    return dot(v.block(), v.vector(0, false), w.vector(0, false));
  } else if ((v.isreal()) && (w.iscomplex())) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v2.vectorC(0, false), w.vectorC(0, false));
  } else if ((w.isreal()) && (v.iscomplex())) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v.vectorC(0, false), w.vectorC(0, false));
  } else {
    return dot(v.block(), v.vectorC(0, false), w.vectorC(0, false));
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double inner(BondList const &bonds, State const &v) try {
  auto w = v;
  apply(bonds, v, w);
  return dot(w, v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

complex innerC(BondList const &bonds, State const &v) try {
  auto w = v;
  apply(bonds, v, w);
  return dotC(w, v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double inner(Bond const &bond, State const &v) try {
  return inner(BondList({bond}), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

complex innerC(Bond const &bond, State const &v) try {
  return innerC(BondList({bond}), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double inner(State const &v, BondList const &bonds, State const &w) try {
  auto bw = zeros_like(w);
  apply(bonds, w, bw);
  return dot(v, bw);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

complex innerC(State const &v, BondList const &bonds, State const &w) try {
  auto bw = zeros_like(w);
  apply(bonds, w, bw);
  return dotC(v, bw);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

double inner(State const &v, Bond const &bond, State const &w) try {
  return inner(v, BondList({bond}), w);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

complex innerC(State const &v, Bond const &bond, State const &w) try {
  return innerC(v, BondList({bond}), w);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

State &operator*=(State &X, complex alpha) try {
  if (X.isreal()) {
    X.make_complex();
    X.matrixC(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return X;
}

State &operator*=(State &X, double alpha) try {
  if (X.isreal()) {
    X.matrix(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return X;
}

State &operator/=(State &X, complex alpha) try {
  if (X.isreal()) {
    X.make_complex();
    X.matrixC(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return X;
}

State &operator/=(State &X, double alpha) try {
  if (X.isreal()) {
    X.matrix(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
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

double dot(block_variant_t const &block, arma::vec const &v,
           arma::vec const &w) try {
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return arma::dot(v, w);
#ifdef XDIAG_USE_MPI
  }
#endif
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

complex dot(block_variant_t const &block, arma::cx_vec const &v,
            arma::cx_vec const &w) try {
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    return cdot_distributed(v, w);
  } else {
#else
  (void)block;
#endif
    return arma::cdot(v, w);
#ifdef XDIAG_USE_MPI
  }
#endif
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

template <typename coeff_t>
double norm(block_variant_t const &block, arma::Col<coeff_t> const &v) try {
  return std::sqrt(xdiag::real(dot(block, v, v)));
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

template double norm(block_variant_t const &, arma::Col<double> const &);
template double norm(block_variant_t const &, arma::Col<complex> const &);

template <typename coeff_t>
double norm1(block_variant_t const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, 1);
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return nrm;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

template double norm1(block_variant_t const &, arma::Col<double> const &);
template double norm1(block_variant_t const &, arma::Col<complex> const &);

template <typename coeff_t>
double norminf(block_variant_t const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, "inf");
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  return nrm;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
  return 0;
}

template double norminf(block_variant_t const &, arma::Col<double> const &);
template double norminf(block_variant_t const &, arma::Col<complex> const &);

} // namespace xdiag
