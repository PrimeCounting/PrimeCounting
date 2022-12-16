#include "mobius_using_newton.h"

#include "../NTT/ntt.h"
#include "../helpers/cell.h"
#include "../helpers/indicators.h"
#include "../helpers/sieve_primes.h"

namespace mobius::details {
size_t get_max_power(prime_t upto, prime_t base) {
  prime_t v = 1;
  size_t ans = 0;
  while (v <= upto / base) {
    v *= base;
    ++ans;
  }
  return ans;
}

std::vector<mint> get_mobius_prime_range(prime_t upto, double lg2_prec,
                                         prime_t min_prime, prime_t max_prime,
                                         size_t vec_sz) {
  size_t max_power =
      std::min(get_max_power(upto, min_prime), max_number_of_factors(upto));
  std::vector<mint> primes_vec;
  {
    // ComputePrimeVector
    primes_vec.resize(vec_sz);
    auto primes = get_primes_by_sieve(max_prime, min_prime);
    add_as_counter(primes_vec, primes, lg2_prec);
    ntt(primes_vec, "Ntt of primes");
  }

  std::vector<mint> mobius;
  {
    // ComputeMobius
    mobius.resize(vec_sz);
    constexpr size_t max_power_available = 1ull << 4;
    ASSERT_FATAL(max_power < max_power_available);

    mint prime_powers[max_power_available];
    mint unique_mults[max_power_available];
    mint inv_mod[max_power_available];

    for (size_t i = 1; i < max_power_available; ++i) {
      inv_mod[i] = mint(i).inverse();
    }

    prime_powers[0] = unique_mults[0] = 1;  // Only 1.

    const std::string title = std::string("Mobius (") +
                              std::to_string(min_prime) + ", " +
                              std::to_string(max_prime) + ")";
    tqdm::Title tq(title);
    for (size_t i = 0; i < vec_sz; ++i) {
      {
        for (size_t power = 1; power <= max_power; ++power)
          // We want the i-th fft coef from the array f' where a prime that
          // was supposed to be in cell c, appears instead in the cell
          // c*power. That is equivalent to: f'(w^i) = f(w^(i*power)).
          prime_powers[power] = primes_vec[i * power % vec_sz];
      }
      {
        unique_mults[1] = prime_powers[1];
        for (size_t power = 2; power <= max_power; ++power) {
          unique_mults[power] = prime_powers[power];
          for (size_t k = power - 1; k > 0; --k) {
            // Apply Newton's identities.
            // The coef of unique_mults[power-1] * prime_power[1] should be 1.
            unique_mults[power] =
                prime_powers[k] * unique_mults[power - k] - unique_mults[power];
          }
          unique_mults[power] *= inv_mod[power];
        }
      }
      {
        mint cur = 0;
        for (int64_t power = max_power; power >= 0; power--) {
          // We want unique_mults[0] * 1;
          cur = unique_mults[power] - cur;
        }
        mobius[i] = cur;
      }
    }
  }
  return mobius;
}

void truncate_and_resize(std::vector<mint>& v, size_t max_cell, size_t new_sz) {
  intt(v, "Truncate INTT");
  ASSERT_FATAL(max_cell < new_sz);
  for (size_t i = max_cell + 1; i < std::min(v.size(), new_sz); ++i) v[i] = 0;
  v.resize(new_sz);
  ntt(v, "Truncate NTT");
}
}  // namespace mobius::details

std::vector<mint> get_mobius_using_newton(prime_t upto, double lg2_prec,
                                          prime_t max_prime) {
  using namespace mobius::details;
  size_t max_cell = get_cell(upto, lg2_prec);
  size_t mobius_sz = ceil_power_of_2(max_cell * 2);
  // Setting mobius_sz to max_cell * 2 also change the thresholds jumps:
  // For 2 this means [p, p^2, p^4, ...] while for 3 it means [p, p^3, p^9, ...]
  // While this does reduce the number of iterations, we pay more per iteration
  // because of the larger vector size (both in the fft and the newton
  // identities).

  constexpr prime_t max_prime_for_naive_conv = 1ll << 10;

  std::vector<prime_t> thresholds({max_prime_for_naive_conv + 1});
  while (thresholds.back() < max_prime + 1) {
    size_t max_power = get_max_power(upto, thresholds.back());
    size_t max_cell_no_overflow = mobius_sz / max_power;
    prime_t prime = get_cell_end(max_cell_no_overflow, lg2_prec);
    thresholds.push_back(std::min(prime, max_prime + 1));
  }

  std::vector<mint> mobius;
  if (thresholds.size() < 2) {
    mobius.resize(mobius_sz);
    mobius[0] = 1;  // {1, 0, 0, 0, ...}
  } else {
    for (size_t i = 1; i < thresholds.size(); ++i) {
      prime_t max_prime_ = thresholds[i] - 1, min_prime_ = thresholds[i - 1];
      size_t max_prime_cell = get_cell(max_prime_, lg2_prec);
      size_t max_power = get_max_power(upto, min_prime_);
      size_t inner_max_cell = max_prime_cell * max_power;
      size_t vec_sz = ceil_power_of_2(inner_max_cell);
      auto cur = get_mobius_prime_range(upto, lg2_prec, min_prime_, max_prime_,
                                        vec_sz);
      if (i == 1) {
        mobius = std::move(cur);  // No need to multiply or truncate.
      } else {
        truncate_and_resize(cur, max_cell, mobius_sz);
        truncate_and_resize(mobius, max_cell, mobius_sz);
        for (size_t i = 0; i < mobius.size(); ++i) mobius[i] *= cur[i];
      }
    }

    intt(mobius, "INTT finalize mobius");
  }
  for (size_t i = max_cell + 1; i < mobius.size(); ++i) mobius[i] = 0;
  {
    // SmallPrimeNaiveConvolution
    auto small_prime_real_bound = std::min(max_prime_for_naive_conv, max_prime);
    auto small_primes = get_primes_by_sieve(small_prime_real_bound);
    for (size_t p_idx : tqdm::title_range<size_t>("SmallPrimeNaiveConvolution",
                                                  small_primes.size())) {
      prime_t p = small_primes[p_idx];
      int64_t cell = get_cell(p, lg2_prec);
      for (int64_t i = max_cell; i >= cell; --i) {
        mobius[i] -= mobius[i - cell];
      }
    }
  }
  return mobius;
}