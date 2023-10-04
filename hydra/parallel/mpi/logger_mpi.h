#pragma once
#ifdef HYDRA_USE_MPI

#include <mpi.h>

#include <iostream>
#include <string>

#define FMT_HEADER_ONLY
#include "extern/fmt/format.h"

namespace hydra {

class LoggerMPI {
public:
  LoggerMPI() : verbosity_(0){};

  void set_verbosity(int verbosity) { verbosity_ = verbosity; }
  int verbosity() { return verbosity_; }

  template <typename... Args>
  void out(const std::string &format, const Args &...args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
      std::cout << fmt::format(format, args...) << "\n";
  }

  template <typename... Args>
  void warn(const std::string &format, const Args &...args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
      std::cout << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  void err(const std::string &format, const Args &...args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
      std::cerr << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  void out(int level, const std::string &format, const Args &...args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ((rank == 0) && (level <= verbosity_))
      std::cout << fmt::format(format, args...) << "\n";
  }

  template <typename... Args>
  void warn(int level, const std::string &format, const Args &...args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ((rank == 0) && (level <= verbosity_))
      std::cout << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  void err(int level, const std::string &format, const Args &...args) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ((rank == 0) && (level <= verbosity_))
      std::cerr << fmt::format(format, args...) << "\n" << std::flush;
  }
  
  template <typename... Args>
  void operator()(const std::string &format, const Args &...args) {
    out(format, args...);
  }
  
  template <typename... Args>
  void operator()(int level, const std::string &format, const Args &...args) {
    out(level, format, args...);
  }

private:
  int verbosity_;
};

inline LoggerMPI LogMPI;
  
} // namespace hydra

#endif
