#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/to_string.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;

int main() try {
  using namespace xdiag::bits;
  using namespace xdiag::combinatorics;

  constexpr bool time_index = false;

  // Subsets
  {
    int n = 16;
    int64_t cnt = 0;
    auto subs = Subsets<uint64_t>(n);
    auto t0 = rightnow();
    for (auto s : subs) {
      if constexpr (time_index)
        if (subs.index(s) != cnt) {
          Log("{} {}", to_string(s, n), subs.index(s));
        }
      ++cnt;
      // Log("{}", to_string(make_bitset(s), n));
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("Subsets: {} {} {} ", cnt, t, t / cnt);
  }

  // Combinations
  {
    int n = 24;
    int k = 12;
    int64_t cnt = 0;
    auto combs = Combinations<uint64_t>(n, k);
    auto t0 = rightnow();
    for (auto s : combs) {
      if constexpr (time_index)

        if (combs.index(s) != cnt) {
          Log("{} {}", to_string(s, n), combs.index(s));
        }
      ++cnt;
      // Log("{}", to_string(make_bitset(s), n));
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("Combinations (short): {} {} {} ", cnt, t, t / cnt);
  }

  // LinTable
  {
    int n = 24;
    int k = 12;
    int64_t cnt = 0;
    auto combs = LinTable<uint64_t>(n, k);
    auto t0 = rightnow();
    for (auto s : combs) {
      if constexpr (time_index)

        if (combs.index(s) != cnt) {
          Log("{} {}", to_string(s, n), combs.index(s));
        }
      ++cnt;
      // Log("{}", to_string(make_bitset(s), n));
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("LinTable: {} {} {} ", cnt, t, t / cnt);
  }

  {
    int n = 80;
    int k = 4;
    int64_t cnt = 0;
    auto combs = Combinations<Bitset<uint64_t, 2>>(n, k);
    auto t0 = rightnow();
    for (auto s : combs) {
      if constexpr (time_index)

        if (combs.index(s) != cnt) {
          Log("{} {}", to_string(s, n), combs.index(s));
        }
      ++cnt;
      // Log("{}", to_string(s, n));
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("Combinations (long): {} {} {} ", cnt, t, t / cnt);
  }

  {
    int64_t n = 12;
    int64_t q = 4;
    int64_t cnt = 0;
    using A = BitArray<uint64_t, 2>;
    auto bms = BoundedMultisets<A>(n, q);
    auto t0 = rightnow();
    for (auto s : bms) {
      if constexpr (time_index)

        if (bms.index(s) != cnt) {
          Log("{} {}", to_string(s, n), bms.index(s));
        }
      ++cnt;
      // Log("{}", to_string(s, n));
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("BoundedMultisets: {} {} {} ", cnt, t, t / cnt);
  }

  {
    int64_t n = 12;
    int64_t total = 12;
    int64_t q = 4;
    int64_t cnt = 0;
    using A = BitArray<uint64_t, 2>;
    auto bps = BoundedPartitions<A>(n, total, q);
    auto t0 = rightnow();
    for (auto s : bps) {
      if constexpr (time_index)
        if (bps.index(s) != cnt) {
          Log("{} {}", to_string(s, n), bps.index(s));
        }
      ++cnt;
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("BoundedPartitions (short): {} {} {} ", cnt, t, t / cnt);
  }


  {
    int64_t n = 50;
    int64_t total = 4;
    int64_t q = 4;
    int64_t cnt = 0;
    using A = BitArray<BitsetStatic2, 2>;
    auto bps = BoundedPartitions<A>(n, total, q);
    auto t0 = rightnow();
    for (auto s : bps) {
      if constexpr (time_index)
        if (bps.index(s) != cnt) {
          Log("{} {}", to_string(s, n), bps.index(s));
        }
      ++cnt;
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("BoundedPartitions (long): {} {} {} ", cnt, t, t / cnt);
  }


  {
    int64_t n = 12;
    int64_t total = 12;
    int64_t q = 4;
    int64_t cnt = 0;
    using A = BitArray<uint64_t, 2>;
    auto bps = SchaeferTable<A>(n, total, q);
    auto t0 = rightnow();
    for (auto s : bps) {
      if constexpr (time_index)
        if (bps.index(s) != cnt) {
          Log("{} {}", to_string(s, n), bps.index(s));
        }
      ++cnt;
    }
    auto t1 = rightnow();
    auto td = duration_cast<microseconds>(t1 - t0).count();
    double t = (double)td;
    Log("SchaeferTable: {} {} {} ", cnt, t, t / cnt);
  }
  
} catch (Error e) {
  error_trace(e);
}
