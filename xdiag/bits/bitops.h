#pragma once

#include <xdiag/bits/extract_deposit.h>
#include <xdiag/bits/gbit.h>
#include <xdiag/bits/popcnt.h>
#include <string>

// #if defined(__BMI__)
// #warning has BMI1
// #endif

// #if defined(__BMI2__)
// #warning has BMI2
// #endif

namespace xdiag::bits {

// bits_to_string
template <typename bit_t>
std::string bits_to_string(bit_t bits, int n, bool reverse = true);

} // namespace xdiag::bits
