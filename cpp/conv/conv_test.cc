#include <gtest/gtest.h>

#include "naive_conv.h"

using vi = std::vector<int>;
TEST(naive_conv, test_naive_conv) {
  auto arr1 = vi({10, 4, 7});
  auto arr2 = vi({5, 3});
  auto expected = vi({50, 50, 47, 21});

  auto res = naive_conv(arr1, arr2);
  EXPECT_EQ(res, expected);
}
