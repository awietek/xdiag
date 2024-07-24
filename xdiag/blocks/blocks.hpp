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
int64_t n_sites(Block const &block);

bool isreal(Block const &block);
bool isdistributed(Block const &block);

std::ostream &operator<<(std::ostream &out, Block const &block);
std::string to_string(Block const &block);

} // namespace xdiag
