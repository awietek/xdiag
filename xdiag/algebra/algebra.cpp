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
  }

  if ((v.isreal()) && (w.isreal())) {
    return dot(v.vector(0, false), w.vector(0, false));
  } else {
    XDIAG_THROW("Unable to compute real dot product of a complex "
                "state, consider using dotC instead.");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
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
  } else if ((v.isreal()) && (!w.isreal())) {
    State v2;
    try {
      v2 = v;
      v2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v2.vectorC(0, false), w.vectorC(0, false));
  } else if ((w.isreal()) && (!v.isreal())) {
    State w2;
    try {
      w2 = w;
      w2.make_complex();
    } catch (...) {
      XDIAG_THROW("Unable to create intermediate complex State");
    }
    return dot(v.block(), v.vectorC(0, false), w2.vectorC(0, false));
  } else {
    return dot(v.block(), v.vectorC(0, false), w.vectorC(0, false));
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

double inner(OpSum const &ops, State const &v) try {
  if (v.isreal() && ops.isreal()) {
    auto w = v;
    apply(ops, v, w);
    return dot(w, v);
  } else {
    XDIAG_THROW("\"inner\" function computing product <psi | O | psi> can only "
                "be called if both the state and the Ops are real. Maybe use "
                "innerC(...) instead.");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

double inner(Op const &op, State const &v) try {
  return inner(OpSum({op}), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

complex innerC(OpSum const &ops, State const &v) try {
  if (v.isreal() && ops.isreal()) {
    auto w = v;
    apply(ops, v, w);
    return (complex)dot(w, v);
  } else if (v.isreal() && !ops.isreal()) {
    auto v2 = v;
    auto w = v2;
    v2.make_complex();
    apply(ops, v2, w);
    return dotC(w, v);
  } else {
    auto w = v;
    apply(ops, v, w);
    return dotC(w, v);
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

complex innerC(Op const &op, State const &v) try {
  return innerC(OpSum({op}), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

// double inner(State const &v, OpSum const &ops, State const &w) try {
//   auto bw = zeros_like(w);
//   apply(ops, w, bw);
//   return dot(v, bw);
// } catch (Error const &error) {
//   XDIAG_RETHROW(error);
// }

// complex innerC(State const &v, OpSum const &ops, State const &w) try {
//   auto bw = zeros_like(w);
//   apply(ops, w, bw);
//   return dotC(v, bw);
// } catch (Error const &error) {
//   XDIAG_RETHROW(error);
// }

// double inner(State const &v, Op const &op, State const &w) try {
//   return inner(v, OpSum({op}), w);
// } catch (Error const &error) {
//   XDIAG_RETHROW(error);
// }

// complex innerC(State const &v, Op const &op, State const &w) try {
//   return innerC(v, OpSum({op}), w);
// } catch (Error const &error) {
//   XDIAG_RETHROW(error);
// }

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
// coeff_t inner(OpSum const &ops, State<coeff_t> const &v) {
//   auto Hv = State<coeff_t>(v.block());
//   apply(ops, v, Hv);
//   return dot(v, Hv);
// }
// template double inner(OpSum const &, State<double> const &);
// template complex inner(OpSum const &, State<complex> const &);

// template <class coeff_t>
// coeff_t inner(Op const &op, State<coeff_t> const &v) {
//   OpSum ops;
//   ops << op;
//   return inner(ops, v);
// }
// template double inner(Op const &, State<double> const &);
// template complex inner(Op const &, State<complex> const &);

// template <class coeff_t>
// coeff_t inner(State<coeff_t> const &w, OpSum const &ops,
//               State<coeff_t> const &v) {
//   auto Hv = State<coeff_t>(w.block());
//   apply(ops, v, Hv);
//   return dot(w, Hv);
// }
// template double inner(State<double> const &, OpSum const &,
//                       State<double> const &);
// template complex inner(State<complex> const &, OpSum const &,
//                        State<complex> const &);

// template <class coeff_t>
// coeff_t inner(State<coeff_t> const &w, Op const &op,
//               State<coeff_t> const &v) {
//   OpSum ops;
//   ops << op;
//   return inner(w, ops, v);
// }
// template double inner(State<double> const &, Op const &,
//                       State<double> const &);
// template complex inner(State<complex> const &, Op const &,
//                        State<complex> const &);

// template <class coeff_t>
// void apply(OpSum const &ops, Block const &block_in,
//            arma::Col<coeff_t> const &vec_in, Block const &block_out,
//            arma::Col<coeff_t> &vec_out) {

//   std::visit(
//       variant::overloaded{
//           [&ops, &vec_in, &vec_out](Spinhalf const &blk_in,
//                                       Spinhalf const &blk_out) {
//             apply(ops, blk_in, vec_in, blk_out, vec_out);
//           },
//           [&ops, &vec_in, &vec_out](tJ const &blk_in, tJ const &blk_out) {
//             apply(ops, blk_in, vec_in, blk_out, vec_out);
//           },
//           [&ops, &vec_in, &vec_out](Electron const &blk_in,
//                                       Electron const &blk_out) {
//             apply(ops, blk_in, vec_in, blk_out, vec_out);
//           },
//           [](auto const &blk_in, auto const &blk_out) {
//             Log.err("Error in apply: Invalid blocks or combination of
//             blocks"); (void)blk_in; (void)blk_out;
//           },
//       },
//       block_in.variant(), block_out.variant());
// }
// template void apply(OpSum const &, Block const &, arma::Col<double> const
// &,
//                     Block const &, arma::Col<double> &);
// template void apply(OpSum const &, Block const &, arma::Col<complex> const
// &,
//                     Block const &, arma::Col<complex> &);

// template <class coeff_t>
// void apply(OpSum const &ops, State<coeff_t> const &state_in,
//            State<coeff_t> &state_out) {
//   apply(ops, state_in.block(), state_in.vector(), state_out.block(),
//         state_out.vector());
// }
// template void apply(OpSum const &, State<double> const &, State<double>
// &); template void apply(OpSum const &, State<complex> const &,
// State<complex> &);

// template <class coeff_t>
// void apply(Op const &op, State<coeff_t> const &state_in,
//            State<coeff_t> &state_out) {
//   OpSum ops;
//   ops << op;
//   apply(ops, state_in, state_out);
// }
// template void apply(Op const &, State<double> const &, State<double> &);
// template void apply(Op const &, State<complex> const &, State<complex> &);

// // template <class coeff_t>
// // void apply(Op const &op, complex coeff, State<coeff_t> const
// &state_in,
// //            State<coeff_t> &state_out) {
// //   OpSum ops;
// //   ops << op;
// //   apply(ops, state_in, state_out);
// // }

// // template <class coeff_t>
// // State<coeff_t> apply(Op const &op, State<coeff_t> const &state_in) {
// //   OpSum ops;
// //   ops << op;
// //   return apply(ops, state_in);
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

double dot(Block const &block, arma::vec const &v,
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
}

complex dot(Block const &block, arma::cx_vec const &v,
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
}

template <typename coeff_t>
double norm(Block const &block, arma::Col<coeff_t> const &v) try {
  return std::sqrt(xdiag::real(dot(block, v, v)));
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template double norm(Block const &, arma::Col<double> const &);
template double norm(Block const &, arma::Col<complex> const &);

template <typename coeff_t>
double norm1(Block const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, 1);
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return nrm;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template double norm1(Block const &, arma::Col<double> const &);
template double norm1(Block const &, arma::Col<complex> const &);

template <typename coeff_t>
double norminf(Block const &block, arma::Col<coeff_t> const &v) try {
  double nrm = arma::norm(v, "inf");
#ifdef XDIAG_USE_MPI
  if (isdistributed(block)) {
    mpi::Allreduce((double *)MPI_IN_PLACE, &nrm, 1, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  return nrm;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template double norminf(Block const &, arma::Col<double> const &);
template double norminf(Block const &, arma::Col<complex> const &);

} // namespace xdiag
