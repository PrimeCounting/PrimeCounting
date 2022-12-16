#pragma once
#include <concepts>
#include <cstddef>
#include <vector>

#include "assertion.h"
#include "double_int.h"
#include "types.h"

template <std::integral T>
constexpr T pow_mod(T a, size_t p, T mod) {
  T res = 1;
  for (; p >= 1; p /= 2, a = wide_mul(a, a) % mod) {
    if (p & 1) {
      res = wide_mul(res, a) % mod;
    }
  }
  return res;
}

prime_t inv_mod(prime_t a, prime_t mod);

template <std::integral T>
constexpr bool naive_is_prime(T v) {
  for (T i = 2; i * i <= v; ++i)
    if (v % i == 0) return false;
  return true;
}

template <std::integral T>
constexpr bool is_prime(T v) {
  if (v <= 3) return true;
  T ord = v - 1;
  size_t num_2 = 0;
  while (ord % 2 == 0) ord /= 2, ++num_2;
  for (size_t test = 2; test <= 7; ++test) {
    T a = test;
    if (a == v) return true;
    a = pow_mod<T>(a, ord, v);
    for (size_t i = 0; a != 1 && i < num_2; ++i) {
      T new_a = wide_mul(a, a) % v;
      if (new_a == 1 && a != v - 1) return false;
      a = new_a;
    }
    if (a != 1) return false;
  }
  return true;
}

template <std::integral T>
std::vector<std::pair<T, size_t>> factor(T t) {
  ASSERT_FATAL(t >= 1);
  std::vector<std::pair<T, size_t>> res;
  for (T v = 2; t > 1; ++v) {
    size_t count = 0;
    while (t % v == 0) {
      t /= v;
      count++;
    }
    if (count) res.emplace_back(v, count);
  }
  return res;
}

constexpr size_t ceil_power_of_2(size_t v) {
  size_t ans = 1;
  while (ans < v) ans *= 2;
  return ans;
}

size_t max_number_of_factors(prime_t v);

template <typename T>
constexpr T pow(T base, size_t power) {
  T res = 1;
  while (power--) res *= base;
  return res;
}
static_assert(pow(2, 5) == 32);
static_assert(pow(3, 2) == 9);
static_assert(pow(100, 0) == 1);

constexpr prime_t pow10(size_t power) { return pow<prime_t>(10, power); }
constexpr prime_t pow2(size_t power) { return pow<prime_t>(2, power); }