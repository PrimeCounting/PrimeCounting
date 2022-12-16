#pragma once
#include <array>
#include <memory>
#include <optional>
#include <vector>

#include "../helpers/assertion.h"
#include "../helpers/double_int.h"
#include "../helpers/indicators.h"
#include "../helpers/types.h"

using Factorization = std::vector<std::pair<prime_t, size_t>>;
/**
 * In order to not save all the factors of the number in a list,
 * we instead do the following:
 *  1. For each value `v` in the range [start, end] we save 2 numbers,
 *     co-prime divisors of `v`, each less then `upto ^ 0.5`, such that
 *     if one would remove all primes that are in the two numbers from
 *     `v` he would remain with a prime number.
 *     We assure this way that each number can be splitted to either
 *     primes or composites no larger than `upto ^ 0.5`.
 *  2. To handle those, we also hold a map from each such number to
 *     one of its factors, thus being able to factor them in
 *     O(number of prime factors).
 *
 * For 1. we sieve on the area, and for each cell we try to add the prime to one
 * of its two divisors (we add if it does not become to large for us to factor).
 */
struct FactorizeBase {
  // `end` is used to calibrate the maximum number to save a factor for.
  FactorizeBase(prime_t end, std::optional<prime_t> max_prime = std::nullopt);

  size_t m_factorize_upto;
  std::vector<prime_t> m_single_factor;

  std::optional<prime_t> m_max_prime;
  std::vector<prime_t> m_primes;
};

struct FactorizedSegment {
  struct SieveCell {
    // Only save upto sqrt of prime_t;
    half_int_t<prime_t> unique_factors_prod[2];
    size_t num_hits;
  };

  FactorizedSegment(prime_t start, std::vector<SieveCell> sieve,
                    std::shared_ptr<FactorizeBase> base);
  void factorize(Factorization& res, prime_t v);
  size_t num_hits_in_sieve(prime_t v);
  const SieveCell& get_value(prime_t v) const;

  prime_t m_start;
  std::vector<SieveCell> m_sieve;
  std::shared_ptr<FactorizeBase> m_base;

  // Sieve the range [start, end) and return the correspondent
  // FactorizedSegment.
  static FactorizedSegment sieve_segment(std::shared_ptr<FactorizeBase>,
                                         prime_t start, prime_t end);
};

template <class ShouldFactorize, class CallBack>
auto handle_range(prime_t start, prime_t end,
                  ShouldFactorize&& should_factorize, CallBack&& call_back,
                  std::optional<std::string> title = std::nullopt,
                  std::optional<prime_t> max_prime = std::nullopt) {
  ASSERT_FATAL(start < end);

  auto base = std::make_shared<FactorizeBase>(end, max_prime);
  // Recommended sieve size:
  constexpr size_t EXTRA_SIEVE_FACTOR = 1ull << 1;
  prime_t segment_size = base->m_primes.back() * EXTRA_SIEVE_FACTOR;
  prime_t num_sieves = (end - start + segment_size - 1) / segment_size;  // ceil

  Factorization factors;
  decltype(call_back(prime_t(1), factors)) ans{};
  std::optional<tqdm::TRange<size_t>> tq;
  if (title.has_value()) {
    tq = tqdm::title_range<size_t>(title.value(), num_sieves);
  }
  // We go in reverse because this makes the progress-bar more indicative.
  for (prime_t i = num_sieves; i > 0; --i) {
    prime_t cur_end = std::min(start + i * segment_size, end);
    prime_t cur_start = start + (i - 1) * segment_size;
    auto fs = FactorizedSegment::sieve_segment(base, cur_start, cur_end);
    for (prime_t v = cur_start; v < cur_end; ++v) {
      if (should_factorize(v, fs.num_hits_in_sieve(v))) {
        fs.factorize(factors, v);
        ans += call_back(v, factors);
      }
    }
    if (tq.has_value()) ++tq.value();
  }
  return ans;
}