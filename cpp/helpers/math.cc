#include "math.h"

#include <array>
#include <limits>

#include "assertion.h"

bool is_prime(prime_t v) {
  if (v <= 3) return true;
  size_t ord = v - 1;
  size_t num_2 = 0;
  while (ord % 2 == 0) ord /= 2, ++num_2;
  constexpr size_t NUM_TESTS = 5;
  for (size_t test = 0; test < NUM_TESTS; ++test) {
    prime_t a;
    do {
      a = rand() % v;
    } while (a == 0);
    a = pow_mod(a, ord, v);
    for (size_t i = 0; a != 1 && i < num_2; ++i) {
      prime_t new_a = wide_mul(a, a) % v;
      if (new_a == 1 && a != v - 1) return false;
      a = new_a;
    }
    if (a != 1) return false;
  }
  return true;
}

namespace {
consteval size_t max_num_factors() {
  prime_t x = std::numeric_limits<prime_t>::max();
  size_t num_primes = 0;
  for (prime_t i = 2; i <= x; ++i) {
    if (naive_is_prime(i)) {
      x /= i;
      ++num_primes;
    }
  }
  return num_primes;
}

consteval auto get_primes_cum_prod() {
  std::array<prime_t, max_num_factors() + 1> res;
  res[0] = 1;
  size_t num_primes = 1;
  for (prime_t i = 2; num_primes < res.size(); ++i) {
    if (naive_is_prime(i)) {
      res[num_primes] = res[num_primes - 1] * i;
      ++num_primes;
    }
  }
  return res;
}
}  // namespace

size_t max_number_of_factors(prime_t v) {
  constexpr auto cum_prod_primes = get_primes_cum_prod();
  ASSERT_FATAL(v > 0);

  for (size_t num_primes = 0; num_primes < cum_prod_primes.size();
       ++num_primes) {
    if (v < cum_prod_primes[num_primes]) {
      return num_primes - 1;
    }
  }
  return cum_prod_primes.size() - 1;
}