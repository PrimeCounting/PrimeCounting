#include "factorize_range.h"

#include <gtest/gtest.h>

#include <iostream>

#include "../helpers/benchmark.h"
#include "../helpers/math.h"
#include "../helpers/types.h"

namespace {
auto as_set(const auto& v) { return std::set(v.begin(), v.end()); }
const auto TRUE_FILTER = [](auto...) { return true; };
}  // namespace

TEST(factorize_range, test_correctness) {
  constexpr prime_t start = 1'000;
  constexpr prime_t end = start * 20;
  auto test = [](prime_t v, const Factorization& factors) -> int {
    EXPECT_EQ(as_set(factors), as_set(factor(v)));
    return 0;
  };
  handle_range(start, end, TRUE_FILTER, test);
}

TEST(factorize_range, test_correctness_with_max_prime) {
  constexpr prime_t start = 1'000;
  constexpr prime_t end = start * 20;
  constexpr prime_t max_prime = 400;
  auto g_max_prime = [&](const auto& v) { return v.first > max_prime; };

  auto test = [&g_max_prime](prime_t v, const Factorization& factors) -> int {
    auto expected = factor(v);
    std::erase_if(expected, g_max_prime);
    EXPECT_EQ(as_set(factors), as_set(expected));
    return 0;
  };
  handle_range(start, end, TRUE_FILTER, test, std::nullopt, max_prime);
}

TEST(benchmark, factorize_base_constructor) {
  constexpr prime_t start = 1ll << 48;
  constexpr prime_t end = start + (1ll << 26);
  auto f = [&]() {
    FactorizeBase fr(end);
    return fr.m_single_factor[50];
  };
  benchmark_code(f, 10);
}

TEST(benchmark, segment_sieve) {
  constexpr prime_t start = 1ll << 48;
  constexpr prime_t end = start + (1ll << 26);
  auto f = [&]() {
    auto base = std::make_shared<FactorizeBase>(end);
    auto fs = FactorizedSegment::sieve_segment(base, start, end);
    size_t ans = 0;
    for (size_t i = start; i < end; ++i) {
      ans += fs.num_hits_in_sieve(i);
    }
    return ans;
  };
  benchmark_code(f, 10);
}

TEST(benchmark, factorize_range_factorize) {
  constexpr prime_t start = 1ll << 44;
  constexpr prime_t num_reps = (1ll << 24);
  constexpr prime_t end = start + num_reps;
  auto base = std::make_shared<FactorizeBase>(end);
  auto fr = FactorizedSegment::sieve_segment(base, start, end);
  auto f = [&]() {
    Factorization factors;
    size_t res = 0;
    for (prime_t i = start; i < end; ++i) {
      fr.factorize(factors, i);
      for (auto v : factors) {
        // To disallow the compiler to optimize the code:
        res += v.first;
      }
    }
    return res;
  };
  benchmark_code(f, 20);
}

TEST(benchmark, factorize_range_all) {
  constexpr prime_t start = 1ll << 44;
  constexpr prime_t num_reps = (1ll << 24);
  constexpr prime_t end = start + num_reps;
  auto f = [&]() {
    auto base = std::make_shared<FactorizeBase>(end);
    return handle_range(start, end, TRUE_FILTER,
                        [&](prime_t, const Factorization& factors) {
                          prime_t res = 0;
                          for (auto v : factors) {
                            // To disallow the compiler to optimize the code:
                            res += v.first;
                          }
                          return res;
                        });
  };
  std::cout << "Start bench" << std::endl;
  benchmark_code(f, 10);
}
