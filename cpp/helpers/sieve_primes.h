#pragma once
#include <optional>
#include <vector>

#include "types.h"

std::vector<prime_t> get_primes_by_sieve(
    prime_t upto, std::optional<prime_t> min_prime = std::nullopt);
