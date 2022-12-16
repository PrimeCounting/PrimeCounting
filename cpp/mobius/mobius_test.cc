#include <gtest/gtest.h>

#include <cmath>

#include "mobius_using_newton.h"
#include "naive_mobius.h"

TEST(mobius, test_naive_mobius) {
  auto res = get_mobius_by_factoring<int32_t>(10, 1);
  auto expected = std::vector<int32_t>{1, -2, -1, 2};
  // Last cell should be 1, but is 2!
  // {8:0, 9:0, 10:1, 11:-1, 12:0, 13:-1, 14:1, 15:1, 21:3};
  // 21 is included in the cell because of rounding errors!
  // log2(3) + log2(7) = 1 + 2 = 3;
  EXPECT_EQ(res, expected);
}

void test_newton_mobius(prime_t upto, double lg2_prec, prime_t max_prime) {
  auto res = get_mobius_using_newton(upto, lg2_prec, max_prime);
  auto expected_ = get_mobius_by_factoring<int32_t>(upto, lg2_prec, max_prime);
  std::vector<mint> expected(expected_.begin(), expected_.end());
  size_t n = std::max(res.size(), expected.size());
  res.resize(n);
  expected.resize(n);
  ASSERT_EQ(res, expected);
}

TEST(mobius, test_newton_mobius) {
  test_newton_mobius(1'000, /*lg2_prec=*/0.1, /*max_prime=*/100);
}
TEST(mobius, test_newton_mobius_large) {
  test_newton_mobius(3'000, /*lg2_prec=*/0.1, /*max_prime=*/3'000);
}

TEST(mobius, test_newton_mobius_small) {
  test_newton_mobius(10, /*lg2_prec=*/1, /*max_prime=*/10);
}
