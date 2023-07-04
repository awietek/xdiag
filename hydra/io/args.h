#pragma once

#include <string>
#include <map>
#include <initializer_list>
#include <utility>

#include <hydra/io/args_handler.h>

namespace hydra {

class Args {
public:
  Args() = default;
  Args(std::initializer_list<std::pair<std::string, io::args_t>> const& args);

  bool defined(std::string key) const;
  io::ArgsHandler operator[](std::string key);

  bool operator==(Args const &other) const;
  bool operator!=(Args const &other) const;

private:
  std::map<std::string, io::args_t> args_;
};

} // namespace hydra
