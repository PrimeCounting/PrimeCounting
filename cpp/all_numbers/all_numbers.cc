#include "all_numbers.h"

#include "../helpers/cell.h"

std::vector<prime_t> get_all_numbers(prime_t upto, double lg2_prec) {
  size_t max_cell = get_cell(upto, lg2_prec);
  std::vector<prime_t> res(max_cell + 1);
  prime_t last_cell_end = 0;
  for (size_t i = 0; i <= max_cell; ++i) {
    prime_t cell_end = get_cell_end(i, lg2_prec);
    res[i] = cell_end - last_cell_end;
    last_cell_end = cell_end;
  }
  return res;
}