#pragma once
#include <cmath>
#include <vector>

#include "types.h"

// Notice that those are the only functions that should use floating points!
// Numeric errors that are caused by this are acceptable.
// The important thing is that this function is monotone.
namespace inexact_functions {
inline size_t get_cell(double v, double lg2_prec) {
  auto res = std::log2(v) / lg2_prec;
  return std::floor(res);
}

inline std::pair<size_t, double> get_cell_and_rem(double v, double lg2_prec) {
  auto res = std::log2(v) / lg2_prec;
  size_t idx = std::floor(res);
  return {idx, res - idx};
}
}  // namespace inexact_functions
using namespace inexact_functions;

prime_t get_cell_end(size_t cell, double lg2_prec);

void add_as_counter(auto& container, const std::vector<prime_t>& primes,
                    double lg2_prec) {
  for (auto prime : primes) {
    size_t cell = get_cell(prime, lg2_prec);
    ++container.at(cell);
  }
}

size_t get_max_value_to_check(size_t max_cell, double lg2_prec);