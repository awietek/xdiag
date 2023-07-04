#include "args_handler.h"

#include <hydra/utils/logger.h>

namespace hydra::io {

ArgsHandler::ArgsHandler(std::string key,
                         std::map<std::string, io::args_t> &args)
    : key_(key), args_(args) {}

template <class data_t> data_t ArgsHandler::as() const {
  if (args_.count(key_)) {
    return std::get<data_t>(args_.at(key_));
  } else {
    Log.err("Error using Args: key \"{}\" not found", key_);
    return data_t();
  }
}
template <class data_t> data_t ArgsHandler::as(data_t const &defaultt) const {
  if (args_.count(key_)) {
    return std::get<data_t>(args_.at(key_));
  } else {
    return defaultt;
  }
}

template <class data_t> void ArgsHandler::operator=(data_t const &data) {
  args_[key_] = data;
}

template bool ArgsHandler::as<bool>() const;
template std::string ArgsHandler::as<std::string>() const;
template int8_t ArgsHandler::as<int8_t>() const;
template int16_t ArgsHandler::as<int16_t>() const;
template int32_t ArgsHandler::as<int32_t>() const;
template int64_t ArgsHandler::as<int64_t>() const;
template uint8_t ArgsHandler::as<uint8_t>() const;
template uint16_t ArgsHandler::as<uint16_t>() const;
template uint32_t ArgsHandler::as<uint32_t>() const;
template uint64_t ArgsHandler::as<uint64_t>() const;
template double ArgsHandler::as<double>() const;
template complex ArgsHandler::as<complex>() const;

template bool ArgsHandler::as(bool const &) const;
template std::string ArgsHandler::as(std::string const &) const;
template int8_t ArgsHandler::as(int8_t const &) const;
template int16_t ArgsHandler::as(int16_t const &) const;
template int32_t ArgsHandler::as(int32_t const &) const;
template int64_t ArgsHandler::as(int64_t const &) const;
template uint8_t ArgsHandler::as(uint8_t const &) const;
template uint16_t ArgsHandler::as(uint16_t const &) const;
template uint32_t ArgsHandler::as(uint32_t const &) const;
template uint64_t ArgsHandler::as(uint64_t const &) const;
template double ArgsHandler::as(double const &) const;
template complex ArgsHandler::as(complex const &) const;

template void ArgsHandler::operator=(bool const &);
template void ArgsHandler::operator=(std::string const &);
template void ArgsHandler::operator=(int8_t const &);
template void ArgsHandler::operator=(int16_t const &);
template void ArgsHandler::operator=(int32_t const &);
template void ArgsHandler::operator=(int64_t const &);
template void ArgsHandler::operator=(uint8_t const &);
template void ArgsHandler::operator=(uint16_t const &);
template void ArgsHandler::operator=(uint32_t const &);
template void ArgsHandler::operator=(uint64_t const &);
template void ArgsHandler::operator=(double const &);
template void ArgsHandler::operator=(complex const &);

} // namespace hydra::io
