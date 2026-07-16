// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <map>
#include <string>
#include <utility>

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_cdagc_string.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_exchange.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_hopdn.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_hopup.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_raise_lower.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_szsz.hpp>
#include <xdiag/math/binomial.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>

using namespace xdiag;

// Apply an elementary operator string (Cup/Cdagup/Cdn/Cdagdn), right-to-left, in
// full electron space with "all ups then all dns" Jordan-Wigner signs (dn ops
// additionally pick up (-1)^Nup for the current ups). Returns valid=false if any
// operator annihilates a missing / creates an existing particle.
struct StringResult {
  bool valid;
  double sign;
  uint32_t ups;
  uint32_t dns;
};
static StringResult
apply_string(std::vector<std::pair<std::string, int64_t>> const &ops,
             uint32_t ups, uint32_t dns, int64_t nsites) {
  double sign = 1.0;
  for (int64_t k = (int64_t)ops.size() - 1; k >= 0; --k) {
    std::string const &t = ops[k].first;
    int64_t s = ops[k].second;
    uint32_t bit = bits::zero<uint32_t>(nsites);
    bits::set(bit, s);
    uint32_t blw = bits::bitmask<uint32_t>(nsites, s);
    if (t == "Cup") {
      if (!bits::get(ups, s)) {
        return {false, 0.0, 0u, 0u};
      }
      if (bits::popcount(ups & blw) & 1) {
        sign = -sign;
      }
      ups ^= bit;
    } else if (t == "Cdagup") {
      if (bits::get(ups, s)) {
        return {false, 0.0, 0u, 0u};
      }
      if (bits::popcount(ups & blw) & 1) {
        sign = -sign;
      }
      ups |= bit;
    } else if (t == "Cdn") {
      if (!bits::get(dns, s)) {
        return {false, 0.0, 0u, 0u};
      }
      if ((bits::popcount(dns & blw) + bits::popcount(ups)) & 1) {
        sign = -sign;
      }
      dns ^= bit;
    } else { // Cdagdn
      if (bits::get(dns, s)) {
        return {false, 0.0, 0u, 0u};
      }
      if ((bits::popcount(dns & blw) + bits::popcount(ups)) & 1) {
        sign = -sign;
      }
      dns |= bit;
    }
  }
  return {true, sign, ups, dns};
}

// Validate the option-B compressed-space dn hop kernel against a brute-force
// full-space dn hop. The kernel computes the Jordan-Wigner sign from the
// COMPRESSED dn string (dnc & fermimask on compressed ranks); the brute force
// computes it from the FULL dn string (dns & between(s1,s2)). Agreement of the
// two matrices is exactly the statement that the compressed computation is
// correct.
TEST_CASE("tj_term_hopdn", "[tj]") {
  using bit_t = uint32_t;
  using enum_t = combinatorics::LinTable<bit_t>;

  for (int64_t nsites = 1; nsites <= 5; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= nsites; ++ndn) {
        enum_t enum_up(nsites, nup);
        enum_t enum_dncs(nsites - nup, ndn);
        basis::BasistJ<enum_t> b(enum_up, enum_dncs);
        int64_t dim = b.size();

        // (ups, full dns) -> linear index, in iteration order (== kernel index)
        std::map<std::pair<bit_t, bit_t>, int64_t> index;
        {
          int64_t idx = 0;
          for (auto [ups, dns] : b) {
            index[{ups, dns}] = idx++;
          }
        }

        for (int64_t s1 = 0; s1 < nsites; ++s1) {
          for (int64_t s2 = 0; s2 < nsites; ++s2) {
            if (s1 == s2) {
              continue;
            }

            // --- kernel matrix ---
            arma::mat mat(dim, dim, arma::fill::zeros);
            auto fill = [&](int64_t idx_in, int64_t idx_out, double val) {
              mat(idx_out, idx_in) += val;
            };
            kernels::tj::term_hopdn<double>(
                Coeff(1.0), Op("Hopdn", std::vector<int64_t>{s1, s2}), b, b,
                fill);

            // --- brute-force reference (full dn string) ---
            arma::mat ref(dim, dim, arma::fill::zeros);
            int64_t l = std::min(s1, s2);
            int64_t u = std::max(s1, s2);
            bit_t flip = bits::zero<bit_t>(nsites);
            bits::set(flip, s1);
            bits::set(flip, s2);
            bit_t between = bits::bitmask<bit_t>(nsites, u - l - 1) << (l + 1);
            int64_t idx_in = 0;
            for (auto [ups, dns] : b) {
              // a dn cannot sit on / hop through an up-occupied site
              if (!bits::get(ups, s1) && !bits::get(ups, s2)) {
                if (bits::popcount(dns & flip) & 1) { // one dn end occupied
                  bit_t dns_out = dns ^ flip;
                  bool fermi = bits::popcount(dns & between) & 1;
                  int64_t idx_out = index.at({ups, dns_out});
                  ref(idx_out, idx_in) += fermi ? 1.0 : -1.0;
                }
              }
              ++idx_in;
            }

            REQUIRE(arma::approx_equal(mat, ref, "absdiff", 1e-12));

            // === Hopup: up hop changes ups -> dn re-mapped (slot surgery) ===
            arma::mat mat_up(dim, dim, arma::fill::zeros);
            auto fill_up = [&](int64_t idx_in, int64_t idx_out, double val) {
              mat_up(idx_out, idx_in) += val;
            };
            kernels::tj::term_hopup<double>(
                Coeff(1.0), Op("Hopup", std::vector<int64_t>{s1, s2}), b, b,
                fill_up);

            // brute force: up hops in full space, dn is an unchanged spectator;
            // the output is in the basis only if the destination is dn-empty.
            arma::mat ref_up(dim, dim, arma::fill::zeros);
            idx_in = 0;
            for (auto [ups, dns] : b) {
              bool up1 = bits::get(ups, s1);
              bool up2 = bits::get(ups, s2);
              if (up1 != up2) {                  // exactly one up among {s1,s2}
                int64_t dest = up1 ? s2 : s1;    // up moves here
                if (!bits::get(dns, dest)) {     // destination free of a dn
                  bit_t ups_out = ups ^ flip;
                  bool fermi = bits::popcount(ups & between) & 1;
                  int64_t idx_out = index.at({ups_out, dns});
                  ref_up(idx_out, idx_in) += fermi ? 1.0 : -1.0;
                }
              }
              ++idx_in;
            }

            REQUIRE(arma::approx_equal(mat_up, ref_up, "absdiff", 1e-12));

            // === SzSz and tJSzSz: diagonal, compressed ranks hoisted ===
            arma::mat mat_sz(dim, dim, arma::fill::zeros);
            arma::mat mat_tjsz(dim, dim, arma::fill::zeros);
            auto fill_sz = [&](int64_t idx_in, int64_t idx_out, double val) {
              mat_sz(idx_out, idx_in) += val;
            };
            auto fill_tjsz = [&](int64_t idx_in, int64_t idx_out, double val) {
              mat_tjsz(idx_out, idx_in) += val;
            };
            kernels::tj::term_szsz<double>(
                Coeff(1.0), Op("SzSz", std::vector<int64_t>{s1, s2}), b, b,
                fill_sz);
            kernels::tj::term_szsz<double>(
                Coeff(1.0), Op("tJSzSz", std::vector<int64_t>{s1, s2}), b, b,
                fill_tjsz);

            arma::mat ref_sz(dim, dim, arma::fill::zeros);
            arma::mat ref_tjsz(dim, dim, arma::fill::zeros);
            idx_in = 0;
            for (auto [ups, dns] : b) {
              int64_t two_sz_i =
                  (int64_t)bits::get(ups, s1) - (int64_t)bits::get(dns, s1);
              int64_t two_sz_j =
                  (int64_t)bits::get(ups, s2) - (int64_t)bits::get(dns, s2);
              ref_sz(idx_in, idx_in) += 0.25 * (double)(two_sz_i * two_sz_j);
              // tJSzSz = SzSz - 1/4 n_i n_j: 0 (parallel), -1/2 (antiparallel)
              int64_t n_i = std::abs(two_sz_i), n_j = std::abs(two_sz_j);
              ref_tjsz(idx_in, idx_in) +=
                  0.25 * (double)(two_sz_i * two_sz_j) - 0.25 * (double)(n_i * n_j);
              ++idx_in;
            }

            REQUIRE(arma::approx_equal(mat_sz, ref_sz, "absdiff", 1e-12));
            REQUIRE(arma::approx_equal(mat_tjsz, ref_tjsz, "absdiff", 1e-12));
          }
        }
      }
    }
  }
}

// Validate the four single creation/annihilation kernels (which change Nup or
// Ndn, so basis_in != basis_out -> rectangular matrices) against full-space
// single-fermion application with the standard "all ups then all dns" Jordan-
// Wigner ordering.
TEST_CASE("tj_raise_lower_kernels", "[tj]") {
  using bit_t = uint32_t;
  using enum_t = combinatorics::LinTable<bit_t>;

  auto make_index = [](basis::BasistJ<enum_t> const &b) {
    std::map<std::pair<bit_t, bit_t>, int64_t> index;
    int64_t idx = 0;
    for (auto [ups, dns] : b) {
      index[{ups, dns}] = idx++;
    }
    return index;
  };

  for (int64_t nsites = 1; nsites <= 5; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= nsites; ++ndn) {
        enum_t eu(nsites, nup), ed(nsites - nup, ndn);
        basis::BasistJ<enum_t> bin(eu, ed);
        auto index_in = make_index(bin);

        for (int64_t s = 0; s < nsites; ++s) {
          bit_t mask = bits::zero<bit_t>(nsites);
          bits::set(mask, s);
          bit_t below = bits::bitmask<bit_t>(nsites, s);

          // op kind -> (nup_out, ndn_out) and validity
          struct Kind {
            const char *name;
            int64_t dnup, dndn;
          };
          for (Kind k : {Kind{"Cup", -1, 0}, Kind{"Cdagup", +1, 0},
                         Kind{"Cdn", 0, -1}, Kind{"Cdagdn", 0, +1}}) {
            int64_t nuo = nup + k.dnup, ndo = ndn + k.dndn;
            if (nuo < 0 || ndo < 0 || nuo + ndo > nsites) {
              continue;
            }
            enum_t euo(nsites, nuo), edo(nsites - nuo, ndo);
            basis::BasistJ<enum_t> bout(euo, edo);
            auto index_out = make_index(bout);
            int64_t dim_in = bin.size(), dim_out = bout.size();

            arma::mat mat(dim_out, dim_in, arma::fill::zeros);
            auto fill = [&](int64_t idx_in, int64_t idx_out, double val) {
              mat(idx_out, idx_in) += val;
            };
            Op op(k.name, s);
            std::string name = k.name;
            if (name == "Cup") {
              kernels::tj::term_cup<double>(Coeff(1.0), op, bin, bout, fill);
            } else if (name == "Cdagup") {
              kernels::tj::term_cdagup<double>(Coeff(1.0), op, bin, bout, fill);
            } else if (name == "Cdn") {
              kernels::tj::term_cdn<double>(Coeff(1.0), op, bin, bout, fill);
            } else {
              kernels::tj::term_cdagdn<double>(Coeff(1.0), op, bin, bout, fill);
            }

            arma::mat ref(dim_out, dim_in, arma::fill::zeros);
            int64_t idx_in = 0;
            for (auto [ups, dns] : bin) {
              bool up_s = bits::get(ups, s);
              bool dn_s = bits::get(dns, s);
              bool neg = bits::popcount(ups) & 1; // (-1)^Nup for dn operators
              if (name == "Cup" && up_s) {
                bool fermi = bits::popcount(ups & below) & 1;
                int64_t io = index_out.at({ups ^ mask, dns});
                ref(io, idx_in) += fermi ? -1.0 : 1.0;
              } else if (name == "Cdagup" && !up_s && !dn_s) {
                bool fermi = bits::popcount(ups & below) & 1;
                int64_t io = index_out.at({ups | mask, dns});
                ref(io, idx_in) += fermi ? -1.0 : 1.0;
              } else if (name == "Cdn" && !up_s && dn_s) {
                bool fermi = (bits::popcount(dns & below) & 1) ^ neg;
                int64_t io = index_out.at({ups, dns ^ mask});
                ref(io, idx_in) += fermi ? -1.0 : 1.0;
              } else if (name == "Cdagdn" && !up_s && !dn_s) {
                bool fermi = (bits::popcount(dns & below) & 1) ^ neg;
                int64_t io = index_out.at({ups, dns | mask});
                ref(io, idx_in) += fermi ? -1.0 : 1.0;
              }
              ++idx_in;
            }

            INFO(std::string(k.name) + " nsites=" + std::to_string(nsites) +
                 " nup=" + std::to_string(nup) + " ndn=" + std::to_string(ndn) +
                 " s=" + std::to_string(s));
            REQUIRE(arma::approx_equal(mat, ref, "absdiff", 1e-12));
          }
        }
      }
    }
  }
}

// Validate the general term_cdagc_string (the cold-path oracle) against an
// independent full-space application of the same normal-ordered operator string:
// each elementary operator is applied in turn with explicit "all ups then all
// dns" Jordan-Wigner signs (dn operators additionally pick up (-1)^Nup for the
// current ups), and any output with double occupancy is projected out.
TEST_CASE("tj_cdagc_string", "[tj]") {
  using bit_t = uint32_t;
  using enum_t = combinatorics::LinTable<bit_t>;
  using OpSeq = std::vector<std::pair<std::string, int64_t>>;

  auto make_index = [](basis::BasistJ<enum_t> const &b) {
    std::map<std::pair<bit_t, bit_t>, int64_t> index;
    int64_t idx = 0;
    for (auto [ups, dns] : b) {
      index[{ups, dns}] = idx++;
    }
    return index;
  };

  for (int64_t nsites = 1; nsites <= 4; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= nsites; ++ndn) {
        enum_t eu(nsites, nup), ed(nsites - nup, ndn);
        basis::BasistJ<enum_t> bin(eu, ed);
        auto index_in = make_index(bin);

        for (int64_t a = 0; a < nsites; ++a) {
          for (int64_t b = 0; b < nsites; ++b) {
            if (a == b) {
              continue;
            }
            // normal-ordered (creation-left, per sector) test strings
            std::vector<OpSeq> monos = {
                {{"Cup", a}},
                {{"Cdagup", a}},
                {{"Cdn", a}},
                {{"Cdagdn", a}},
                {{"Cdagup", a}, {"Cup", b}},                       // up hop
                {{"Cdagdn", a}, {"Cdn", b}},                       // dn hop
                {{"Cdagup", a}, {"Cup", b}, {"Cdagdn", b}, {"Cdn", a}}, // mixed
            };

            for (OpSeq const &ops : monos) {
              int64_t dnup = 0, dndn = 0;
              std::vector<Op> opv;
              for (auto const &[t, s] : ops) {
                opv.push_back(Op(t, s));
                if (t == "Cdagup") {
                  ++dnup;
                } else if (t == "Cup") {
                  --dnup;
                } else if (t == "Cdagdn") {
                  ++dndn;
                } else {
                  --dndn;
                }
              }
              int64_t nuo = nup + dnup, ndo = ndn + dndn;
              if (nuo < 0 || ndo < 0 || nuo + ndo > nsites) {
                continue;
              }
              enum_t euo(nsites, nuo), edo(nsites - nuo, ndo);
              basis::BasistJ<enum_t> bout(euo, edo);
              auto index_out = make_index(bout);
              int64_t dim_in = bin.size(), dim_out = bout.size();

              arma::mat mat(dim_out, dim_in, arma::fill::zeros);
              auto fill = [&](int64_t i, int64_t o, double v) { mat(o, i) += v; };
              kernels::tj::term_cdagc_string<double>(Coeff(1.0), Monomial(opv),
                                                      bin, bout, fill);

              // full-space reference: apply ops right-to-left with JW signs
              arma::mat ref(dim_out, dim_in, arma::fill::zeros);
              int64_t idx_in = 0;
              for (auto [ups, dns] : bin) {
                bit_t u = ups, d = dns;
                double sign = 1.0;
                bool valid = true;
                for (int64_t k = (int64_t)ops.size() - 1; k >= 0 && valid; --k) {
                  std::string const &t = ops[k].first;
                  int64_t s = ops[k].second;
                  bit_t bit = bits::zero<bit_t>(nsites);
                  bits::set(bit, s);
                  bit_t blw = bits::bitmask<bit_t>(nsites, s);
                  if (t == "Cup") {
                    if (!bits::get(u, s)) { valid = false; break; }
                    if (bits::popcount(u & blw) & 1) { sign = -sign; }
                    u ^= bit;
                  } else if (t == "Cdagup") {
                    if (bits::get(u, s)) { valid = false; break; }
                    if (bits::popcount(u & blw) & 1) { sign = -sign; }
                    u |= bit;
                  } else if (t == "Cdn") {
                    if (!bits::get(d, s)) { valid = false; break; }
                    int64_t f = bits::popcount(d & blw) + bits::popcount(u);
                    if (f & 1) { sign = -sign; }
                    d ^= bit;
                  } else { // Cdagdn
                    if (bits::get(d, s)) { valid = false; break; }
                    int64_t f = bits::popcount(d & blw) + bits::popcount(u);
                    if (f & 1) { sign = -sign; }
                    d |= bit;
                  }
                }
                if (valid && ((u & d) == 0u)) { // project out double occupancy
                  ref(index_out.at({u, d}), idx_in) += sign;
                }
                ++idx_in;
              }

              REQUIRE(arma::approx_equal(mat, ref, "absdiff", 1e-12));
            }
          }
        }
      }
    }
  }
}

// Validate the fast term_exchange against the full-space Exchange operator
// S+{i} S-{j} + S-{i} S+{j} (each a four-fermion string applied with explicit
// Jordan-Wigner signs, double occupancy projected out). Exchange conserves Nup
// and Ndn, so the matrix is square on the same block.
TEST_CASE("tj_exchange", "[tj]") {
  using bit_t = uint32_t;
  using enum_t = combinatorics::LinTable<bit_t>;
  using OpSeq = std::vector<std::pair<std::string, int64_t>>;

  for (int64_t nsites = 1; nsites <= 5; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= nsites; ++ndn) {
        enum_t eu(nsites, nup), ed(nsites - nup, ndn);
        basis::BasistJ<enum_t> b(eu, ed);
        int64_t dim = b.size();

        std::map<std::pair<bit_t, bit_t>, int64_t> index;
        {
          int64_t idx = 0;
          for (auto [ups, dns] : b) {
            index[{ups, dns}] = idx++;
          }
        }

        for (int64_t i = 0; i < nsites; ++i) {
          for (int64_t j = 0; j < nsites; ++j) {
            if (i == j) {
              continue;
            }

            arma::mat mat(dim, dim, arma::fill::zeros);
            auto fill = [&](int64_t in, int64_t o, double v) { mat(o, in) += v; };
            kernels::tj::term_exchange<double>(
                Coeff(1.0), Op("Exchange", std::vector<int64_t>{i, j}), b, b,
                fill);

            // reference: S+{i}S-{j} + S-{i}S+{j} in full space
            OpSeq spm = {{"Cdagup", i}, {"Cdn", i}, {"Cdagdn", j}, {"Cup", j}};
            OpSeq smp = {{"Cdagdn", i}, {"Cup", i}, {"Cdagup", j}, {"Cdn", j}};
            // Exchange = 1/2 (S+S- + S-S+)
            arma::mat ref(dim, dim, arma::fill::zeros);
            int64_t idx_in = 0;
            for (auto [ups, dns] : b) {
              for (OpSeq const &str : {spm, smp}) {
                StringResult r = apply_string(str, ups, dns, nsites);
                if (r.valid && ((r.ups & r.dns) == 0u)) {
                  ref(index.at({r.ups, r.dns}), idx_in) += 0.5 * r.sign;
                }
              }
              ++idx_in;
            }

            INFO("nsites=" + std::to_string(nsites) + " nup=" +
                 std::to_string(nup) + " ndn=" + std::to_string(ndn) +
                 " i=" + std::to_string(i) + " j=" + std::to_string(j));
            REQUIRE(arma::approx_equal(mat, ref, "absdiff", 1e-12));
          }
        }
      }
    }
  }
}

// The same kernels on the NON-number-conserving (no-np) basis: ups runs over all
// subsets and the dn fiber varies with popcount(ups). The whole space (dim
// 3^nsites) is one block, so basis_in == basis_out for every operator (including
// the number-changing raise/lower ops). This exercises the no-np seam of BasistJ
// (popcount-indexed fibers, prefix-sum offsets) through every kernel category.
TEST_CASE("tj_nonp_kernels", "[tj]") {
  using bit_t = uint32_t;
  using enum_t = combinatorics::Subsets<bit_t>;

  for (int64_t nsites = 1; nsites <= 4; ++nsites) {
    enum_t eu(nsites);
    basis::BasistJ<enum_t> b(eu); // no-np ctor
    int64_t dim = b.size();

    int64_t expect = 1;
    for (int64_t k = 0; k < nsites; ++k) {
      expect *= 3;
    }
    REQUIRE(dim == expect); // 3^nsites

    std::map<std::pair<bit_t, bit_t>, int64_t> index;
    {
      int64_t idx = 0;
      for (auto [ups, dns] : b) {
        index[{ups, dns}] = idx++;
      }
    }

    auto string_ref = [&](std::vector<std::pair<std::string, int64_t>> const
                              &ops) {
      arma::mat r(dim, dim, arma::fill::zeros);
      int64_t in = 0;
      for (auto [ups, dns] : b) {
        StringResult res = apply_string(ops, ups, dns, nsites);
        if (res.valid && ((res.ups & res.dns) == 0u)) {
          r(index.at({res.ups, res.dns}), in) += res.sign;
        }
        ++in;
      }
      return r;
    };

    // hops, szsz, exchange, and the general string, per bond
    for (int64_t i = 0; i < nsites; ++i) {
      for (int64_t j = 0; j < nsites; ++j) {
        if (i == j) {
          continue;
        }
        int64_t l = std::min(i, j), u = std::max(i, j);
        bit_t flip = bits::zero<bit_t>(nsites);
        bits::set(flip, i);
        bits::set(flip, j);
        bit_t between = bits::bitmask<bit_t>(nsites, u - l - 1) << (l + 1);

        // --- Hopdn ---
        arma::mat mhd(dim, dim, arma::fill::zeros);
        auto fhd = [&](int64_t in, int64_t o, double v) { mhd(o, in) += v; };
        kernels::tj::term_hopdn<double>(
            Coeff(1.0), Op("Hopdn", std::vector<int64_t>{i, j}), b, b, fhd);
        arma::mat rhd(dim, dim, arma::fill::zeros);
        // --- Hopup ---
        arma::mat mhu(dim, dim, arma::fill::zeros);
        auto fhu = [&](int64_t in, int64_t o, double v) { mhu(o, in) += v; };
        kernels::tj::term_hopup<double>(
            Coeff(1.0), Op("Hopup", std::vector<int64_t>{i, j}), b, b, fhu);
        arma::mat rhu(dim, dim, arma::fill::zeros);
        // --- tJSzSz ---
        arma::mat msz(dim, dim, arma::fill::zeros);
        auto fsz = [&](int64_t in, int64_t o, double v) { msz(o, in) += v; };
        kernels::tj::term_szsz<double>(
            Coeff(1.0), Op("SzSz", std::vector<int64_t>{i, j}), b, b, fsz);
        arma::mat rsz(dim, dim, arma::fill::zeros);

        int64_t in = 0;
        for (auto [ups, dns] : b) {
          // hopdn ref
          if (!bits::get(ups, i) && !bits::get(ups, j)) {
            if (bits::popcount(dns & flip) & 1) {
              bit_t dns_out = dns ^ flip;
              bool f = bits::popcount(dns & between) & 1;
              rhd(index.at({ups, dns_out}), in) += f ? 1.0 : -1.0;
            }
          }
          // hopup ref
          if (bits::get(ups, i) != bits::get(ups, j)) {
            int64_t dest = bits::get(ups, i) ? j : i;
            if (!bits::get(dns, dest)) {
              bool f = bits::popcount(ups & between) & 1;
              rhu(index.at({ups ^ flip, dns}), in) += f ? 1.0 : -1.0;
            }
          }
          // szsz ref
          int64_t szi = (int64_t)bits::get(ups, i) - (int64_t)bits::get(dns, i);
          int64_t szj = (int64_t)bits::get(ups, j) - (int64_t)bits::get(dns, j);
          rsz(in, in) += 0.25 * (double)(szi * szj);
          ++in;
        }
        REQUIRE(arma::approx_equal(mhd, rhd, "absdiff", 1e-12));
        REQUIRE(arma::approx_equal(mhu, rhu, "absdiff", 1e-12));
        REQUIRE(arma::approx_equal(msz, rsz, "absdiff", 1e-12));

        // --- Exchange ---
        arma::mat mex(dim, dim, arma::fill::zeros);
        auto fex = [&](int64_t in2, int64_t o, double v) { mex(o, in2) += v; };
        kernels::tj::term_exchange<double>(
            Coeff(1.0), Op("Exchange", std::vector<int64_t>{i, j}), b, b, fex);
        arma::mat rex = 0.5 * (string_ref({{"Cdagup", i}, {"Cdn", i},
                                           {"Cdagdn", j}, {"Cup", j}}) +
                               string_ref({{"Cdagdn", i}, {"Cup", i},
                                           {"Cdagup", j}, {"Cdn", j}}));
        REQUIRE(arma::approx_equal(mex, rex, "absdiff", 1e-12));

        // --- general Cdag/C string (mixed 4-op) ---
        arma::mat mcs(dim, dim, arma::fill::zeros);
        auto fcs = [&](int64_t in2, int64_t o, double v) { mcs(o, in2) += v; };
        std::vector<Op> opv = {Op("Cdagup", i), Op("Cup", j), Op("Cdagdn", j),
                               Op("Cdn", i)};
        kernels::tj::term_cdagc_string<double>(Coeff(1.0), Monomial(opv), b, b,
                                                fcs);
        arma::mat rcs = string_ref(
            {{"Cdagup", i}, {"Cup", j}, {"Cdagdn", j}, {"Cdn", i}});
        REQUIRE(arma::approx_equal(mcs, rcs, "absdiff", 1e-12));
      }
    }

    // raise/lower single operators (Nup/Ndn change but stay in the same block)
    for (int64_t s = 0; s < nsites; ++s) {
      arma::mat m0(dim, dim, arma::fill::zeros), m1(dim, dim, arma::fill::zeros),
          m2(dim, dim, arma::fill::zeros), m3(dim, dim, arma::fill::zeros);
      auto f0 = [&](int64_t in, int64_t o, double v) { m0(o, in) += v; };
      auto f1 = [&](int64_t in, int64_t o, double v) { m1(o, in) += v; };
      auto f2 = [&](int64_t in, int64_t o, double v) { m2(o, in) += v; };
      auto f3 = [&](int64_t in, int64_t o, double v) { m3(o, in) += v; };
      kernels::tj::term_cup<double>(Coeff(1.0), Op("Cup", s), b, b, f0);
      kernels::tj::term_cdagup<double>(Coeff(1.0), Op("Cdagup", s), b, b, f1);
      kernels::tj::term_cdn<double>(Coeff(1.0), Op("Cdn", s), b, b, f2);
      kernels::tj::term_cdagdn<double>(Coeff(1.0), Op("Cdagdn", s), b, b, f3);
      REQUIRE(arma::approx_equal(m0, string_ref({{"Cup", s}}), "absdiff", 1e-12));
      REQUIRE(arma::approx_equal(m1, string_ref({{"Cdagup", s}}), "absdiff", 1e-12));
      REQUIRE(arma::approx_equal(m2, string_ref({{"Cdn", s}}), "absdiff", 1e-12));
      REQUIRE(arma::approx_equal(m3, string_ref({{"Cdagdn", s}}), "absdiff", 1e-12));
    }
  }
}
