#include "say_hello.hpp"

#include <iostream>
#include <sstream>

#include <xdiag/config.hpp>

#ifndef XDIAG_DISABLE_COLOR
#include <xdiag/extern/fmt/color.hpp>
#endif

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
#ifdef XDIAG_DISABLE_COLOR
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
#else
  // // standard
  // auto c1 = fg(fmt::color::blue);
  // auto c2 = fg(fmt::color::orange);
  // auto c3 = fg(fmt::color::red);
  // auto c4 = fg(fmt::color::magenta);
  // auto c5 = fg(fmt::color::green);

  // spring
  auto c1 = fg(fmt::rgb(0xAFDDD5));
  auto c2 = fg(fmt::rgb(0xFFA700));
  auto c3 = fg(fmt::rgb(0xFFCCCD));
  auto c4 = fg(fmt::rgb(0xF56093));
  auto c5 = fg(fmt::rgb(0x64864A));
  
  std::vector<std::string> lines(9);

  lines[0] = std::string(" _____     __   ______     ") +
             fmt::format(c3, R"(.-./`))") + std::string("    ____      .-") +
             fmt::format(c5, "_") + std::string("'''-.");

  lines[1] = std::string(R"( \   )") + fmt::format(c1, "_") +
             std::string(R"(\   /  / |    )") + fmt::format(c2, "_") +
             std::string(R"( `''. )") + fmt::format(c3, R"(\ .-.'))") +
             std::string(R"( .'  __ `.  ')") + fmt::format(c5, "_( )_") +
             std::string(R"(   \)");

  lines[2] = std::string(R"( )") + fmt::format(c1, ".-./ )") +
             std::string(R"(. /  '  | )") + fmt::format(c2, "_ | ) _") +
             std::string(R"(  \)") + fmt::format(c3, R"(/ `-' \)") +
             std::string(R"(/   '  \  \|)") + fmt::format(c5, "(_ o _)") +
             std::string(R"(|  ')");

  lines[3] = std::string(R"( )") + fmt::format(c1, R"(\ '_ .'))") +
             std::string(R"( .'   |)") + fmt::format(c2, "( ''_'  )") +
             std::string(R"( | )") + fmt::format(c3, R"(`-'`"`)") +
             std::string(R"(|___|  /  |. )") + fmt::format(c5, "(_,_)") +
             std::string(R"(/___|)");

  lines[4] = fmt::format(c1, R"((_ (_) _))") + std::string(R"( '    | )") +
             fmt::format(c2, ". (_) `.") +
             std::string(R"( | .---.    _.-`   ||  |  .-----.)");

  lines[5] = std::string(R"(  )") + fmt::format(c1, R"(/    \)") +
             std::string(R"(   \   |)") + fmt::format(c2, "(_    ._)") +
             std::string(R"( ' |   | .'   )") + fmt::format(c4, "_") +
             std::string(R"(    |'  \  '-   .')");

  lines[6] = std::string(R"(  )") + fmt::format(c1, R"(`-'`-')") +
             std::string(R"(    \  |  )") + fmt::format(c2, R"((_.\.')") +
             std::string(R"( /  |   | |  )") + fmt::format(c4, "_( )_") +
             std::string(R"(  | \  `-'`   |)");

  lines[7] = std::string(R"( /  /   \    \ |       .'   |   | \ )") +
             fmt::format(c4, "(_ o _)") + std::string(R"( /  \        /)");

  lines[8] = std::string(R"('--'     '----''-----'`     '---'  '.)") +
             fmt::format(c4, "(_,_)") + std::string(R"(.'    `'-...-')");

  //   std::string hello =
  //       R"(
  //  _____     __   ______     .-./`)    ____      .-_'''-.
  //  \   _\   /  / |    _ `''. \ .-.') .'  __ `.  '_( )_   \
//  .-./ ). /  '  | _ | ) _  \/ `-' \/   '  \  \|(_ o _)|  '
  //  \ '_ .') .'   |( ''_'  ) | `-'`"`|___|  /  |. (_,_)/___|
  // (_ (_) _) '    | . (_) `. | .---.    _.-`   ||  |  .-----.
  //   /    \   \   |(_    ._) ' |   | .'   _    |'  \  '-   .'
  //   `-'`-'    \  |  (_.\.' /  |   | |  _( )_  | \  `-'`   |
  //  /  /   \    \ |       .'   |   | \ (_ o _) /  \        /
  // '--'     '----''-----'`     '---'  '.(_,_).'    `'-...-'
  // )";
  std::string hello = "";
  for (auto line : lines) {
    hello += line + std::string("\n");
  }
#endif
  std::cout << hello << "\n" << version_string() << "\n";
}

void print_version() { std::cout << version_string() << "\n"; }
} // namespace xdiag
