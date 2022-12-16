#pragma once
#include <cmath>
#include <optional>
#include <vector>

#include "../helpers/cell.h"
#include "../helpers/types.h"
#include "naive_mobius.h"

template <std::integral T>
std::vector<T> get_mobius_by_factoring(
    prime_t upto, double lg2_prec,
    std::optional<prime_t> max_prime = std::nullopt) {
  size_t max_cell = get_cell(upto, lg2_prec);
  prime_t max_value_to_check = get_max_value_to_check(max_cell, lg2_prec);

  std::vector<T> res(max_cell + 1);
  for (prime_t i = 1; i < max_value_to_check; ++i) {
    size_t cell = 0;
    int mobius = 1;
    for (auto [p, c] : factor(i)) {
      if (c > 1 || (max_prime.has_value() && p > max_prime)) {
        mobius = 0;
        break;
      }
      mobius *= -1;
      cell += get_cell(p, lg2_prec);
    }
    if (cell <= max_cell && mobius != 0) {
      res[cell] += mobius;
    }
  }
  return res;
}
