#include "logger.hpp"

namespace xdiag {

void set_verbosity(int64_t level) { xdiag::Log.set_verbosity(level); }

} // namespace xdiag
