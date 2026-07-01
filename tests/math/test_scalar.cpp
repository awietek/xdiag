// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/math/scalar.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("test_scalar", "[complex]") {
  Log("Testing Scalar");

  // --- Construction ---
  {
    Scalar d(1.5);
    REQUIRE(d.isreal());
    REQUIRE(d.real() == 1.5);
    REQUIRE(d.imag() == 0.0);

    Scalar c(complex(1.0, 2.0));
    REQUIRE(!c.isreal());
    REQUIRE(c.real() == 1.0);
    REQUIRE(c.imag() == 2.0);

    // Integer constructors resolve ambiguity
    Scalar i(2);
    REQUIRE(i.isreal());
    REQUIRE(i.real() == 2.0);

    Scalar i64((int64_t)3);
    REQUIRE(i64.isreal());
    REQUIRE(i64.real() == 3.0);
  }

  // --- Type query: is<T> and as<T> ---
  {
    Scalar d(3.14);
    REQUIRE(d.is<double>());
    REQUIRE(!d.is<complex>());
    REQUIRE(d.as<double>() == Approx(3.14));
    REQUIRE(d.as<complex>() == complex(3.14, 0.0));

    Scalar c(complex(1.0, -1.0));
    REQUIRE(!c.is<double>());
    REQUIRE(c.is<complex>());
    REQUIRE(c.as<complex>() == complex(1.0, -1.0));
    // REQUIRE_THROWS(c.as<double>());
  }

  // --- abs and conj ---
  {
    Scalar d(3.0);
    REQUIRE(d.abs() == Approx(3.0));
    REQUIRE(d.conj().real() == Approx(3.0));
    REQUIRE(d.conj().imag() == Approx(0.0));

    Scalar c(complex(3.0, 4.0));
    REQUIRE(c.abs() == Approx(5.0));
    REQUIRE(c.conj().real() == Approx(3.0));
    REQUIRE(c.conj().imag() == Approx(-4.0));
  }

  // --- Arithmetic: stays real when possible ---
  {
    Scalar a(2.0), b(3.0);
    REQUIRE((a + b).isreal());
    REQUIRE((a + b).real() == Approx(5.0));
    REQUIRE((a - b).real() == Approx(-1.0));
    REQUIRE((a * b).real() == Approx(6.0));
    REQUIRE((a / b).real() == Approx(2.0 / 3.0));
    REQUIRE((-a).real() == Approx(-2.0));
  }

  // --- Arithmetic: widens to complex on mixed input ---
  {
    Scalar r(2.0);
    Scalar c(complex(1.0, 1.0));
    Scalar sum = r + c;
    REQUIRE(!sum.isreal());
    REQUIRE(sum.real() == Approx(3.0));
    REQUIRE(sum.imag() == Approx(1.0));

    Scalar prod = r * c;
    REQUIRE(!prod.isreal());
    REQUIRE(prod.real() == Approx(2.0));
    REQUIRE(prod.imag() == Approx(2.0));
  }

  // --- Compound assignment ---
  {
    Scalar a(1.0);
    a += Scalar(2.0);
    REQUIRE(a.real() == Approx(3.0));
    a *= Scalar(complex(0.0, 1.0));
    REQUIRE(!a.isreal());
    REQUIRE(a.imag() == Approx(3.0));
  }

  // --- zero ---
  {
    Scalar zr = zero(Scalar(5.0));
    REQUIRE(zr.isreal());
    REQUIRE(zr.real() == Approx(0.0));

    Scalar zc = zero(Scalar(complex(1.0, 1.0)));
    REQUIRE(!zc.isreal());
    REQUIRE(zc.abs() == Approx(0.0));
  }

  // --- to_real ---
  {
    Scalar r(2.5);
    REQUIRE(r.to_real().isreal());
    REQUIRE(r.to_real().real() == Approx(2.5));

    Scalar pure_real_complex(complex(3.0, 0.0));
    REQUIRE(pure_real_complex.to_real(1e-10).isreal());
    REQUIRE(pure_real_complex.to_real(1e-10).real() == Approx(3.0));

    Scalar c(complex(1.0, 1.0));
    // REQUIRE_THROWS(c.to_real());
  }

  // --- isapprox ---
  {
    Scalar a(1.0), b(1.0 + 1e-14);
    REQUIRE(isapprox(a, b));

    Scalar c1(complex(1.0, 1.0)), c2(complex(1.0 + 1e-14, 1.0 - 1e-14));
    REQUIRE(isapprox(c1, c2));

    REQUIRE(!isapprox(Scalar(1.0), Scalar(2.0)));
  }

  // --- equality ---
  {
    REQUIRE(Scalar(1.0) == Scalar(1.0));
    REQUIRE(Scalar(1.0) != Scalar(2.0));
    REQUIRE(Scalar(complex(1.0, 2.0)) == Scalar(complex(1.0, 2.0)));
    // double(1.0) != complex(1.0, 0.0): different variant types
    REQUIRE(Scalar(1.0) != Scalar(complex(1.0, 0.0)));
  }
}
