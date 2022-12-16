#pragma once
#include <vector>

#include "../helpers/assertion.h"

template <std::integral T>
std::vector<T> naive_conv(const std::vector<T> &left,
                          const std::vector<T> &right) {
  size_t n1 = left.size(), n2 = right.size();
  ASSERT_FATAL(n1 > 0 and n2 > 0);
  size_t n3 = n1 + n2 - 1;
  std::vector<T> res(n3);
  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < n2; ++j) {
      res[i + j] += left[i] * right[j];
    }
  }
  return res;
}