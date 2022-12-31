#include "all_numbers.h"

#include <gtest/gtest.h>

#include "../helpers/mod_int.h"
#include "../helpers/types.h"
#include "naive_all_numbers.h"

using mi = ModInt32<1'000'000'007>;
using vmi = std::vector<mi>;

TEST(all_numbers, test_naive_all_numbers) {
  auto expected = vmi({1, 2, 4, 8});

  auto res = get_all_numbers_by_for_loop<mi>(10, 1.);
  EXPECT_EQ(res, expected);
}

void test_all_numbers(prime_t upto, double lg2_prec) {
  auto res = get_all_numbers(upto, lg2_prec);
  auto expected = get_all_numbers_by_for_loop<prime_t>(upto, lg2_prec);
  EXPECT_EQ(res, expected);
}

TEST(all_numbers, test_all_numbers) { test_all_numbers(1'000'000, 0.001); }
TEST(all_numbers, test_all_numbers_small) { test_all_numbers(20, 0.5); }
