#include <vector>

#include "../helpers/cell.h"
#include "../helpers/types.h"

template <std::integral T>
std::vector<T> get_all_numbers_by_for_loop(prime_t upto, double lg2_prec) {
  size_t max_cell = get_cell(upto, lg2_prec);
  std::vector<T> res(max_cell + 1);
  for (prime_t i = 1; get_cell(i, lg2_prec) <= max_cell; ++i) {
    ++res[get_cell(i, lg2_prec)];
  }
  return res;
}