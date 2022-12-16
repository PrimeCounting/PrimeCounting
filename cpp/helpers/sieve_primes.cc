#include "sieve_primes.h"

std::vector<prime_t> get_primes_by_sieve(prime_t upto,
                                         std::optional<prime_t> min_prime) {
  std::vector<bool> sieve(upto + 1, true);
  std::vector<prime_t> primes;
  for (prime_t i = 2; i <= upto; ++i) {
    if (sieve[i]) {
      if (min_prime <= i) primes.push_back(i);
      for (prime_t j = i * i; j <= upto; j += i) {
        sieve[j] = false;
      }
    }
  }
  return primes;
}