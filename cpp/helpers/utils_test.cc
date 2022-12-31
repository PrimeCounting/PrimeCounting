#include <gtest/gtest.h>

#include <utility>
#include <vector>

#include "math.h"

TEST(utils_test, is_prime) {
  EXPECT_TRUE(is_prime(2));
  EXPECT_TRUE(is_prime(3));
  EXPECT_TRUE(is_prime(5));
  EXPECT_TRUE(is_prime(101));
  EXPECT_TRUE(is_prime(1'000'000'007));
  EXPECT_FALSE(is_prime(2 * 3));
  EXPECT_FALSE(is_prime(101 * 5));
  EXPECT_FALSE(is_prime(5 * 5));
  for (size_t i = 2; i < 100'000; ++i) {
    ASSERT_EQ(is_prime(i), naive_is_prime(i)) << i;
  }
}

TEST(utils_test, factorize) {
  using R = std::vector<std::pair<int32_t, size_t>>;
  EXPECT_EQ(factor(30), R({{2, 1}, {3, 1}, {5, 1}}));
  EXPECT_EQ(factor(1000), R({{2, 3}, {5, 3}}));
  EXPECT_EQ(factor(13 * 11), R({{11, 1}, {13, 1}}));
  EXPECT_EQ(factor(11 * 11), R({{11, 2}}));
  EXPECT_EQ(factor(1), R());
  EXPECT_ANY_THROW(factor(0));
  EXPECT_ANY_THROW(factor(-1));
  EXPECT_ANY_THROW(factor(-100));
}
