#include "cell.h"

#include "assertion.h"
#include "math.h"

prime_t get_cell_end(size_t cell, double lg2_prec) {
  // We can probably use even less than 10bit range.
  constexpr double machine_prec =
      1. / pow2(std::numeric_limits<double>::digits - 10);
  // This is the cell end approximately.
  double approx_log = (cell + 1) * lg2_prec;
  double delta = approx_log * machine_prec;
  double approx_number = std::exp2(approx_log);
  prime_t l = std::floor(approx_number - approx_number * delta);
  prime_t r = std::ceil(approx_number + approx_number * delta);
  if (l == 0) l = 1;
  ASSERT_FATAL(get_cell(l, lg2_prec) <= cell);
  ASSERT_FATAL(cell < get_cell(r, lg2_prec));
  --r;

  // If this part become a bottle neck, we can improve it
  // by doing `exp2` using long double.
  while (l < r) {
    prime_t m = (l + r + 1) / 2;
    if (get_cell(m, lg2_prec) <= cell) {
      l = m;
    } else {
      r = m - 1;
    }
  }
  return l;
}

size_t get_max_value_to_check(size_t max_cell, double lg2_prec) {
  size_t max_cell_to_check = max_cell;
  while (max_cell_to_check <= max_cell + max_number_of_factors(get_cell_end(
                                             max_cell_to_check, lg2_prec))) {
    ++max_cell_to_check;
  }
  return get_cell_end(max_cell_to_check, lg2_prec);
}