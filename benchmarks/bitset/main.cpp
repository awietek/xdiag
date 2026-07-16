#include <xdiag/utils/timing.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>

using namespace xdiag;

int main() try {
  using combinatorics::Combinations;
  using bits::Bitset;
  
  int n = 34;
  int k = 6;
  int64_t cnt1 = 0;
  tic();
  for (auto s : Combinations<uint64_t>(n, k)) {
    ++cnt1;
  }
  toc();

  int64_t cnt2 = 0;
  tic();
  for (auto s : Combinations<Bitset<uint16_t, 4>>(n, k)) {
    ++cnt2;
  }
  toc();
  Log("{} {}", cnt1, cnt2);

} catch (Error e) {
  error_trace(e);
}
