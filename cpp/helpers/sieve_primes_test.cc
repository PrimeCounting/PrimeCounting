#include "sieve_primes.h"

#include <gtest/gtest.h>

#include "math.h"
#include "types.h"

TEST(primes_classic, naive_test) {
  constexpr prime_t upto = 1'000'000;
  auto primes = get_primes_by_sieve(upto);
  std::vector<prime_t> primes_;
  for (prime_t i = 2; i <= upto; ++i) {
    if (is_prime(i)) {
      primes_.push_back(i);
    }
  }
  EXPECT_EQ(primes, primes_);
}
