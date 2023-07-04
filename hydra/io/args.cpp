#include "args.h"

namespace hydra {

Args::Args(
    std::initializer_list<std::pair<std::string, io::args_t>> const &args) {
  for (auto const &[key, val] : args) {
    args_[key] = val;
  }
}
bool Args::defined(std::string key) const { return args_.count(key); }
io::ArgsHandler Args::operator[](std::string key) {
  return io::ArgsHandler(key, args_);
}

bool Args::operator==(Args const &other) const { return args_ == other.args_; }
bool Args::operator!=(Args const &other) const { return !operator==(other); }

} // namespace hydra
