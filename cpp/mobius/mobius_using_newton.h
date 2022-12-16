#pragma once

#include "../NTT/ntt.h"
#include "../helpers/types.h"

namespace mobius::details {
size_t get_max_power(prime_t upto, prime_t base);
}  // namespace mobius::details

std::vector<mint> get_mobius_using_newton(prime_t upto, double lg2_prec,
                                          prime_t max_prime);
