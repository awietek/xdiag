#include "args.h"

namespace hydra {
Arg::Arg(std::string key, std::string val) : key_(key), val_(val) {}
Arg::Arg(std::string key, const char *val)
    : key_(key), val_(std::string(val)) {}
Arg::Arg(std::string key, bool val) : key_(key), val_(val) {}
Arg::Arg(std::string key, int8_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, int16_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, int32_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, int64_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, uint8_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, uint16_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, uint32_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, uint64_t val) : key_(key), val_(val) {}
Arg::Arg(std::string key, double val) : key_(key), val_(val) {}
Arg::Arg(std::string key, complex val) : key_(key), val_(val) {}

std::string Arg::key() const { return key_; }
io::args_t Arg::val() const { return val_; }

Args::Args(std::initializer_list<Arg> const &args) {
  for (auto const & arg : args) {
    std::string key = arg.key();
    auto val = arg.val();
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
