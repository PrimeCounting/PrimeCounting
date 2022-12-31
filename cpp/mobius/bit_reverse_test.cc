#include "bit_reverse.h"

#include <gtest/gtest.h>

#include <numeric>
#include <vector>

TEST(bit_reverse, test_reverse_array) {
  std::vector<int> v(1 << 16);
  std::iota(v.begin(), v.end(), 0);
  auto v1 = v;
  auto v2 = v;
  auto expected = v1;
  bit_reverse_naive(expected);
  bit_reverse_impl_small(v1);
  inplace_bit_reverse(v2);
  EXPECT_EQ(v1, expected);
  EXPECT_EQ(v2, expected);
}