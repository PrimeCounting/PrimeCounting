#include "ntt.h"

#include <bits/stdc++.h>
#include <gtest/gtest.h>

#include "../conv/naive_conv.h"
#include "../helpers/benchmark.h"

using NT = std::vector<mint>;
TEST(NTT, test_conv) {
  NT arr1({10, 4, 7, 0, 0, 0, 0, 0});
  NT arr2({5, 3, 0, 0, 0, 0, 0, 0});
  ntt(arr1), ntt(arr2);
  NT res = arr1;
  for (size_t i = 0; i < res.size(); ++i) res[i] *= arr2[i];
  intt(res);

  {
    NT expected({50, 50, 47, 21, 0, 0, 0, 0});
    EXPECT_EQ(res, expected);
  }
  for (mint &i : arr1) i *= i;
  intt(arr1);
  {
    NT expected({100, 80, 16 + 140, 2 * 4 * 7, 7 * 7, 0, 0, 0});
    EXPECT_EQ(arr1, expected);
  }
  {
    constexpr size_t n = 1ull << 10;
    NT large(n);
    std::iota(large.begin(), large.begin() + (n / 4), mint(0));
    auto expected = naive_conv(large, large);
    expected.resize(n);
    auto res = large;
    ntt(res);
    for (auto &i : res) i = i * i;
    intt(res);
    EXPECT_EQ(res, expected);
  }
}

TEST(NTT, test_montgomery) {
  mint x = mint::get_mod() - 50;
  EXPECT_EQ(x + 50, 0);
  EXPECT_EQ(x + 51, 1);
  EXPECT_EQ(x - x, 0);
  EXPECT_EQ(x + (-50), -100);
  EXPECT_EQ(x - x - 1, -1);
  EXPECT_EQ(x / x, 1);
}

TEST(benchmark, ntt_conv) {
  NT a(1ull << 25);
  std::iota(a.begin(), a.end(), 1);
  benchmark_code(
      [&]() {
        ntt(a);
        return std::accumulate(a.begin(), a.end(), mint(0));
      },
      15);
}