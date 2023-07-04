#pragma once

#include <string>
#include <map>
#include <initializer_list>
#include <utility>

#include <hydra/io/args_handler.h>

namespace hydra {

class Arg {
public:
  Arg() = default;
  Arg(std::string key, std::string val);
  Arg(std::string key, const char *val);
  Arg(std::string key, bool val);
  Arg(std::string key, int8_t val);
  Arg(std::string key, int16_t val);
  Arg(std::string key, int32_t val);
  Arg(std::string key, int64_t val);
  Arg(std::string key, uint8_t val);
  Arg(std::string key, uint16_t val);
  Arg(std::string key, uint32_t val);
  Arg(std::string key, uint64_t val);
  Arg(std::string key, double val);
  Arg(std::string key, complex val);

  std::string key() const;
  io::args_t val() const;
private:
  std::string key_;
  io::args_t val_;
};
  
class Args {
public:
  Args() = default;
  Args(std::initializer_list<Arg> const& args);

  bool defined(std::string key) const;
  io::ArgsHandler operator[](std::string key);

  bool operator==(Args const &other) const;
  bool operator!=(Args const &other) const;

private:
  std::map<std::string, io::args_t> args_;
};

} // namespace hydra
