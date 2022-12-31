#include "mod_int.h"

#include <gtest/gtest.h>

#include <memory>

#include "benchmark.h"

using mi = ModInt32<1'000'000'007>;

TEST(mod_int, operators) {
  {
    mi x;
    ASSERT_EQ(x, 0);
  }
  {
    mi x = mi::get_mod() + 50;
    ASSERT_EQ(x, 50);
  }
  {
    mi x = mi::get_mod() - 50;
    ASSERT_EQ(x, -50);
    ASSERT_EQ(x.get(), mi::get_mod() - 50);
    ASSERT_EQ(x + 50, 0u);
    ASSERT_EQ(x + 51, 1);
    ASSERT_EQ(x - x, 0);
    ASSERT_EQ(x + (-50), -100);
    ASSERT_EQ(x - x - 1, -1);
    ASSERT_EQ(x / x, 1);
  }
}

TEST(benchmark, mod_int) {
  benchmark_code([]() {
    mi x;
    constexpr size_t n = 100'000'000;
    for (size_t i = 0; i < n; ++i) {
      mi v = rand();
      for (size_t j = 0; j < 10; ++j, ++i) {
        x = (x * v);
        ++x;
      }
    }
    return x.get();
  });
}