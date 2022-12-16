#pragma once
#include <cmath>

#include "../NTT/ntt.h"
#include "../helpers/types.h"
#include "error_correction.h"

mint count_primes_with_errors(prime_t upto, double lg2_prec,
                              prime_t max_prime_to_use);

prime_t count_primes(prime_t upto, double lg2_prec, prime_t max_prime_to_use);

inline prime_t count_primes(prime_t upto) {
  double lg2_prec = 1. / std::sqrt(upto);
  prime_t max_prime_to_use = std::ceil(std::sqrt(upto));
  return count_primes(upto, lg2_prec, max_prime_to_use);
}