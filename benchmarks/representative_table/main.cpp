#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/tables/representative_table.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;

int main() try {
  using bits::Bitset;
  using combinatorics::Combinations;
  using namespace symmetries;

  int n = 36;
  int k = 18;
  auto combs = Combinations<uint64_t>(n, k);
  auto irrep = cyclic_group_irrep(n, 0); // trivial: all chi=1
  auto sp = SitePermutation(irrep.group());
  tic();
  auto table = RepresentativeTable<Combinations<uint64_t>>(combs, sp, irrep);
  toc();

} catch (Error e) {
  error_trace(e);
}
