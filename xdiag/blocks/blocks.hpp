#pragma once
#include <variant>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/common.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#endif

namespace xdiag {

#ifdef XDIAG_USE_MPI
using Block =
    std::variant<Spinhalf, tJ, Electron, SpinhalfDistributed, tJDistributed>;
#else
using Block = std::variant<Spinhalf, tJ, Electron>;
#endif

int64_t dim(Block const &block);
int64_t size(Block const &block);
int64_t nsites(Block const &block);

bool isreal(Block const &block);

constexpr bool isdistributed(Block const &block) {
  return std::visit(
      overload{
          [&](Spinhalf const &) -> bool { return false; },
          [&](tJ const &) -> bool { return false; },
          [&](Electron const &) -> bool { return false; },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &) -> bool { return true; },
          [&](tJDistributed const &) -> bool { return true; },
#endif
          [&](auto &&) -> bool { return false; },
      },
      block);
}
  
std::ostream &operator<<(std::ostream &out, Block const &block);
std::string to_string(Block const &block);

} // namespace xdiag
