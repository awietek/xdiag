// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <sstream>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("state", "[states]") try {
  using namespace xdiag;
  using namespace arma;

  auto block = Spinhalf(4, 2);
  auto psir = State(block);
  REQUIRE(isreal(psir));
  // XDIAG_SHOW(psir);
  // XDIAG_SHOW(psir.matrix());
  REQUIRE(norm(psir.matrix()) < 1e-12);

  auto psic = State(block, false);
  REQUIRE(!isreal(psic));
  // XDIAG_SHOW(psic);
  // XDIAG_SHOW(psic.matrixC());
  REQUIRE(norm(psic.matrixC()) < 1e-12);

  arma::vec vr = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  auto psivr = State(block, vr);
  REQUIRE(isreal(psivr));
  // XDIAG_SHOW(psivr.vector());
  REQUIRE(norm(psivr.vector() - vr) < 1e-12);

  arma::vec vi = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  auto vc = arma::cx_vec(vr, vi);
  auto psivc = State(block, vc);
  // XDIAG_SHOW(psivc.vectorC());
  REQUIRE(!isreal(psivc));
  REQUIRE(norm(psivc.vectorC() - vc) < 1e-12);
  REQUIRE(norm(psivc.real().vector() - vr) < 1e-12);
  REQUIRE(norm(psivc.imag().vector() - vi) < 1e-12);
    
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state matrix constructors", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6

  // --- real multi-column matrix constructor + accessors ---
  {
    arma::mat M(6, 3, arma::fill::randu);
    auto psi = State(block, M);
    REQUIRE(isreal(psi));
    REQUIRE(psi.ncols() == 3);
    REQUIRE(psi.nrows() == 6);
    REQUIRE(arma::norm(psi.matrix() - M) < 1e-12);

    // col(n) extracts a single-column State
    auto c1 = psi.col(1);
    REQUIRE(c1.ncols() == 1);
    REQUIRE(arma::norm(c1.vector() - M.col(1)) < 1e-12);
  }

  // --- complex multi-column matrix constructor ---
  {
    arma::cx_mat M(6, 2, arma::fill::randu);
    auto psi = State(block, M);
    REQUIRE(!isreal(psi));
    REQUIRE(psi.ncols() == 2);
    REQUIRE(arma::norm(psi.matrixC() - M) < 1e-12);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state constructor and accessor errors", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6

  // --- constructor size mismatches throw ---
  REQUIRE_THROWS_AS(State(block, arma::vec(5, arma::fill::zeros)), Error);
  REQUIRE_THROWS_AS(State(block, arma::cx_vec(5, arma::fill::zeros)), Error);
  REQUIRE_THROWS_AS(State(block, arma::mat(5, 2, arma::fill::zeros)), Error);
  REQUIRE_THROWS_AS(State(block, arma::cx_mat(5, 2, arma::fill::zeros)), Error);

  // --- real state: complex accessors throw ---
  {
    auto psi = State(block, true);
    REQUIRE_THROWS_AS(psi.vectorC(), Error);
    REQUIRE_THROWS_AS(psi.matrixC(), Error);
    // out-of-range / negative column indices throw
    REQUIRE_THROWS_AS(psi.vector(5), Error);  // only 1 column
    REQUIRE_THROWS_AS(psi.vector(-1), Error);
  }

  // --- complex state: real accessors throw ---
  {
    auto psi = State(block, false);
    REQUIRE_THROWS_AS(psi.vector(), Error);
    REQUIRE_THROWS_AS(psi.matrix(), Error);
    REQUIRE_THROWS_AS(psi.vectorC(-1), Error);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state accessors and free functions", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6
  arma::mat M(6, 3, arma::fill::randu);
  arma::cx_mat MC(6, 2, arma::fill::randu);
  auto psir = State(block, M);
  auto psic = State(block, MC);

  // --- sizes ---
  REQUIRE(psir.nsites() == 4);
  REQUIRE(nsites(psir) == 4);
  REQUIRE(psir.dim() == 6);
  REQUIRE(dim(psir) == 6);
  REQUIRE(psir.size() == 18);
  REQUIRE(size(psir) == 18);
  REQUIRE(nrows(psir) == 6);
  REQUIRE(ncols(psir) == 3);
  REQUIRE(size(psic) == 12);
  REQUIRE(ncols(psic) == 2);
  REQUIRE(isvalid(psir));
  REQUIRE(!isvalid(State()));

  // --- free-function accessors agree with member accessors ---
  REQUIRE(arma::norm(matrix(psir) - M) < 1e-12);
  REQUIRE(arma::norm(matrixC(psic) - MC) < 1e-12);
  REQUIRE(arma::norm(vector(psir, 2) - M.col(2)) < 1e-12);
  REQUIRE(arma::norm(vectorC(psic, 1) - MC.col(1)) < 1e-12);
  REQUIRE(arma::norm(col(psir, 1).vector() - M.col(1)) < 1e-12);
  REQUIRE(arma::norm(col(psic, 0).vectorC() - MC.col(0)) < 1e-12);
  REQUIRE_THROWS_AS(col(psir, 3), Error);
  REQUIRE_THROWS_AS(col(psic, 2), Error);

  // --- non-copying views alias the storage ---
  {
    auto psi = State(block, M);
    arma::mat view = psi.matrix(false);
    view(0, 0) = 42.0;
    REQUIRE(psi.matrix()(0, 0) == 42.0);
    arma::vec vview = psi.vector(1, false);
    vview(2) = -7.0;
    REQUIRE(psi.matrix()(2, 1) == -7.0);
  }
  {
    auto psi = State(block, MC);
    arma::cx_mat view = psi.matrixC(false);
    view(0, 0) = complex(1.0, 2.0);
    REQUIRE(psi.matrixC()(0, 0) == complex(1.0, 2.0));
  }

  // --- real / imag parts, member and free ---
  {
    REQUIRE(arma::norm(real(psic).matrix() - arma::real(MC)) < 1e-12);
    REQUIRE(arma::norm(imag(psic).matrix() - arma::imag(MC)) < 1e-12);
    REQUIRE(isreal(real(psic)));
    REQUIRE(isreal(imag(psic)));
    REQUIRE(real(psic).ncols() == 2);
    REQUIRE(imag(psic).ncols() == 2);

    // real state: real() is a copy, imag() is zero
    REQUIRE(arma::norm(real(psir).matrix() - M) < 1e-12);
    REQUIRE(arma::norm(imag(psir).matrix()) < 1e-12);
    REQUIRE(imag(psir).ncols() == 3);
  }

  // --- make_complex, member and free ---
  {
    auto psi = State(block, M);
    make_complex(psi);
    REQUIRE(!isreal(psi));
    REQUIRE(arma::norm(arma::real(psi.matrixC()) - M) < 1e-12);
    REQUIRE(arma::norm(arma::imag(psi.matrixC())) < 1e-12);
    // idempotent on a complex state
    make_complex(psi);
    REQUIRE(arma::norm(arma::real(psi.matrixC()) - M) < 1e-12);
  }

  // --- isapprox ---
  {
    auto a = State(block, M);
    arma::mat Mp = M + 1e-14;
    auto b = State(block, Mp);
    REQUIRE(isapprox(a, b));
    arma::mat M2 = M * 2.0;
    REQUIRE(!isapprox(a, State(block, M2)));
    REQUIRE(!isapprox(a, psic)); // real vs complex
    REQUIRE(!isapprox(a, State(Spinhalf(4, 1), true)));
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state developer constructors and pointers", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6
  arma::mat M(6, 2, arma::fill::randu);
  arma::cx_mat MC(6, 2, arma::fill::randu);

  // --- raw-pointer constructors copy the data ---
  {
    auto psi = State(block, M.memptr(), 2);
    REQUIRE(isreal(psi));
    REQUIRE(psi.ncols() == 2);
    REQUIRE(arma::norm(psi.matrix() - M) < 1e-12);

    auto psic = State(block, MC.memptr(), 2);
    REQUIRE(!isreal(psic));
    REQUIRE(arma::norm(psic.matrixC() - MC) < 1e-12);

    // stride 2 reads every other double, i.e. the real parts of MC
    double const *ptr = reinterpret_cast<double const *>(MC.memptr());
    auto re = State(block, ptr, 2, 2);
    REQUIRE(arma::norm(re.matrix() - arma::real(MC)) < 1e-12);
    auto im = State(block, ptr + 1, 2, 2);
    REQUIRE(arma::norm(im.matrix() - arma::imag(MC)) < 1e-12);
    REQUIRE_THROWS_AS(State(block, M.memptr(), 2, 3), Error);
  }

  // --- memptr / colptr expose the column-major storage ---
  {
    auto psi = State(block, M);
    REQUIRE(psi.memptr() == psi.colptr(0));
    REQUIRE(psi.colptr(1) == psi.memptr() + 6);
    REQUIRE(psi.colptr(1)[3] == M(3, 1));
    psi.colptr(1)[3] = 5.0;
    REQUIRE(psi.matrix()(3, 1) == 5.0);
    REQUIRE_THROWS_AS(psi.colptr(2), Error);
    REQUIRE_THROWS_AS(psi.colptr(-1), Error);

    auto psic = State(block, MC);
    REQUIRE(psic.memptrC() == psic.colptrC(0));
    REQUIRE(psic.colptrC(1) == psic.memptrC() + 6);
    REQUIRE(psic.colptrC(1)[2] == MC(2, 1));
    REQUIRE_THROWS_AS(psic.colptrC(2), Error);
    REQUIRE_THROWS_AS(psic.colptrC(-1), Error);
  }

  // --- invalid column count ---
  REQUIRE_THROWS_AS(State(block, true, 0), Error);
  REQUIRE_THROWS_AS(State(block, false, -1), Error);
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state arithmetic", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6
  arma::mat M(6, 2, arma::fill::randu);
  arma::mat N(6, 2, arma::fill::randu);
  arma::cx_mat MC(6, 2, arma::fill::randu);
  arma::cx_mat NC(6, 2, arma::fill::randu);
  complex z(0.5, -1.5);

  // --- scaling by a real scalar ---
  {
    auto a = State(block, M);
    REQUIRE(arma::norm((a * 2.0).matrix() - 2.0 * M) < 1e-12);
    REQUIRE(arma::norm((3.0 * a).matrix() - 3.0 * M) < 1e-12);
    REQUIRE(arma::norm((a / 4.0).matrix() - M / 4.0) < 1e-12);
    a *= 2.0;
    REQUIRE(arma::norm(a.matrix() - 2.0 * M) < 1e-12);
    a /= 2.0;
    REQUIRE(arma::norm(a.matrix() - M) < 1e-12);

    auto c = State(block, MC);
    REQUIRE(arma::norm((c * 2.0).matrixC() - 2.0 * MC) < 1e-12);
    REQUIRE(arma::norm((c / 2.0).matrixC() - MC / 2.0) < 1e-12);
    c *= 2.0;
    c /= 2.0;
    REQUIRE(arma::norm(c.matrixC() - MC) < 1e-12);
  }

  // --- scaling by a complex scalar promotes real states ---
  {
    auto a = State(block, M);
    State az = a * z;
    REQUIRE(!isreal(az));
    REQUIRE(isreal(a)); // operand untouched
    arma::cx_mat expected = arma::conv_to<arma::cx_mat>::from(M) * z;
    REQUIRE(arma::norm(az.matrixC() - expected) < 1e-12);
    REQUIRE(arma::norm((z * a).matrixC() - expected) < 1e-12);
    REQUIRE(arma::norm((a / z).matrixC() -
                       arma::conv_to<arma::cx_mat>::from(M) / z) < 1e-12);

    a *= z;
    REQUIRE(!isreal(a));
    REQUIRE(arma::norm(a.matrixC() - expected) < 1e-12);
    a /= z;
    REQUIRE(arma::norm(arma::real(a.matrixC()) - M) < 1e-12);

    auto c = State(block, MC);
    REQUIRE(arma::norm((c * z).matrixC() - MC * z) < 1e-12);
    REQUIRE(arma::norm((c / z).matrixC() - MC / z) < 1e-12);
    c *= z;
    REQUIRE(arma::norm(c.matrixC() - MC * z) < 1e-12);
    c /= z;
    REQUIRE(arma::norm(c.matrixC() - MC) < 1e-12);
  }

  // --- addition and subtraction in all real/complex combinations ---
  {
    auto a = State(block, M);
    auto b = State(block, N);
    auto c = State(block, MC);
    auto d = State(block, NC);
    arma::cx_mat Mc = arma::conv_to<arma::cx_mat>::from(M);
    arma::cx_mat Nc = arma::conv_to<arma::cx_mat>::from(N);

    REQUIRE(isreal(a + b));
    REQUIRE(arma::norm((a + b).matrix() - (M + N)) < 1e-12);
    REQUIRE(arma::norm((a - b).matrix() - (M - N)) < 1e-12);
    REQUIRE(arma::norm((-a).matrix() + M) < 1e-12);

    REQUIRE(!isreal(a + c));
    REQUIRE(arma::norm((a + c).matrixC() - (Mc + MC)) < 1e-12);
    REQUIRE(arma::norm((c + a).matrixC() - (Mc + MC)) < 1e-12);
    REQUIRE(arma::norm((a - c).matrixC() - (Mc - MC)) < 1e-12);
    REQUIRE(arma::norm((c - a).matrixC() - (MC - Mc)) < 1e-12);
    REQUIRE(arma::norm((c + d).matrixC() - (MC + NC)) < 1e-12);
    REQUIRE(arma::norm((c - d).matrixC() - (MC - NC)) < 1e-12);
    REQUIRE(arma::norm((-c).matrixC() + MC) < 1e-12);

    // compound forms, including promotion of the left operand
    auto x = State(block, M);
    x += b;
    REQUIRE(arma::norm(x.matrix() - (M + N)) < 1e-12);
    x -= b;
    REQUIRE(arma::norm(x.matrix() - M) < 1e-12);
    x += c;
    REQUIRE(!isreal(x));
    REQUIRE(arma::norm(x.matrixC() - (Mc + MC)) < 1e-12);
    x -= c;
    REQUIRE(arma::norm(x.matrixC() - Mc) < 1e-12);
    x -= a;
    REQUIRE(arma::norm(x.matrixC()) < 1e-12);

    auto y = State(block, MC);
    y += a;
    REQUIRE(arma::norm(y.matrixC() - (MC + Mc)) < 1e-12);
    y -= a;
    REQUIRE(arma::norm(y.matrixC() - MC) < 1e-12);
  }

  // --- mismatched operands throw ---
  {
    auto a = State(block, M);
    auto one_col = State(block, arma::vec(6, arma::fill::randu));
    REQUIRE_THROWS_AS(a + one_col, Error);
    REQUIRE_THROWS_AS(a - one_col, Error);
    auto other = State(Spinhalf(4, 1), true, 2);
    REQUIRE_THROWS_AS(a + other, Error);
    REQUIRE_THROWS_AS(a -= other, Error);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state printing", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2);
  {
    auto psi = State(block, true, 3);
    std::string s = to_string(psi);
    REQUIRE(s.find("REAL State") != std::string::npos);
    REQUIRE(s.find("ncols: 3") != std::string::npos);
    std::ostringstream oss;
    oss << psi;
    REQUIRE(oss.str() == s);
  }
  {
    auto psi = State(block, false);
    REQUIRE(to_string(psi).find("COMPLEX State") != std::string::npos);
  }
  {
    REQUIRE(to_string(State()).find("INVALID State") != std::string::npos);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state on every block type", "[states]") try {
  using namespace xdiag;

  std::vector<Block> blocks = {Spinhalf(3, 1), Boson(3, 2, 2), Fermion(3, 1),
                               Electron(3, 1, 1), tJ(3, 1, 1)};
  for (Block const &block : blocks) {
    int64_t D = size(block);
    arma::cx_mat MC(D, 2, arma::fill::randu);
    auto psic = State(block, MC);
    REQUIRE(nsites(psic) == 3);
    REQUIRE(dim(psic) == D);
    REQUIRE(arma::norm(real(psic).matrix() - arma::real(MC)) < 1e-12);
    REQUIRE(arma::norm(imag(psic).matrix() - arma::imag(MC)) < 1e-12);
    REQUIRE(arma::norm(col(psic, 1).vectorC() - MC.col(1)) < 1e-12);

    arma::mat M(D, 2, arma::fill::randu);
    auto psir = State(block, M);
    REQUIRE(arma::norm(col(psir, 1).vector() - M.col(1)) < 1e-12);
    REQUIRE(arma::norm(imag(psir).matrix()) < 1e-12);
    REQUIRE(to_string(psir).find("REAL State") != std::string::npos);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("zero_state and zero", "[states]") try {
  using namespace xdiag;

  std::vector<Block> blocks = {Spinhalf(4, 2), Boson(3, 2, 2), Fermion(4, 2),
                               Electron(3, 1, 1), tJ(3, 1, 1)};
  for (Block const &block : blocks) {
    int64_t D = size(block);

    State zr = zero_state(block);
    REQUIRE(isvalid(zr));
    REQUIRE(isreal(zr));
    REQUIRE(zr.ncols() == 1);
    REQUIRE(zr.nrows() == D);
    REQUIRE(arma::norm(zr.vector()) == 0.0);

    State zc = zero_state(block, false, 3);
    REQUIRE(!isreal(zc));
    REQUIRE(zc.ncols() == 3);
    REQUIRE(arma::norm(zc.matrixC()) == 0.0);

    // zero() clears an existing state in place
    State r = random_state(block, true, 2);
    REQUIRE(arma::norm(r.matrix()) > 0.0);
    zero(r);
    REQUIRE(arma::norm(r.matrix()) == 0.0);
    REQUIRE(r.ncols() == 2);

    State c = random_state(block, false, 2);
    REQUIRE(arma::norm(c.matrixC()) > 0.0);
    zero(c);
    REQUIRE(arma::norm(c.matrixC()) == 0.0);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
