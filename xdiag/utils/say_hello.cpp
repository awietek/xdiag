#include "say_hello.hpp"

#include <iostream>
#include <sstream>

#include <xdiag/config.hpp>

namespace xdiag {

std::string version_string() {
  std::stringstream sstr;
  sstr << "XDiag Version: " << XDIAG_VERSION << "\n";
  sstr << "Git hash     : " << git_hash() << "\n";
  sstr << "Directory    : " << XDIAG_DIRECTORY << "\n";
  sstr << "Hostname     : " << XDIAG_HOSTNAME << "\n";
  sstr << "Compiled by  : " << XDIAG_COMPILEDBY << "\n";
  return sstr.str();
}

void say_hello() {
  std::string hello =
      R"(
 _____     __   ______     .-./`)    ____      .-_'''-.
 \   _\   /  / |    _ `''. \ .-.') .'  __ `.  '_( )_   \
 .-./ ). /  '  | _ | ) _  \/ `-' \/   '  \  \|(_ o _)|  '
 \ '_ .') .'   |( ''_'  ) | `-'`"`|___|  /  |. (_,_)/___|
(_ (_) _) '    | . (_) `. | .---.    _.-`   ||  |  .-----.
  /    \   \   |(_    ._) ' |   | .'   _    |'  \  '-   .'
  `-'`-'    \  |  (_.\.' /  |   | |  _( )_  | \  `-'`   |
 /  /   \    \ |       .'   |   | \ (_ o _) /  \        /
'--'     '----''-----'`     '---'  '.(_,_).'    `'-...-'
)";
  std::cout << hello << "\n" << version_string() << "\n";
}

void print_version() { std::cout << version_string() << "\n"; }
} // namespace xdiag
