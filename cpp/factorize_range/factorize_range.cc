#include "factorize_range.h"

#include <cmath>
#include <iterator>
#include <memory>

#include "../helpers/assertion.h"
#include "../helpers/double_int.h"
#include "../helpers/types.h"

namespace {
std::vector<prime_t> get_single_factor(prime_t upto) {
  std::vector<prime_t> res(upto + 1, 1);
  for (auto i = 2; i <= upto; ++i) {
    if (res[i] == 1) {
      for (auto j = i; j <= upto; j += i) {
        res[j] = i;
      }
    }
  }
  return res;
}
}  // namespace

FactorizeBase::FactorizeBase(prime_t end, std::optional<prime_t> max_prime)
    : m_max_prime(max_prime), m_primes() {
  ASSERT_FATAL(!max_prime.has_value() or max_prime.value() <= end);
  // Note that factorize_upto might be larger than max_prime (and that's ok).
  m_factorize_upto = std::ceil(std::sqrt(end));
  m_single_factor = get_single_factor(m_factorize_upto);

  prime_t max_prime_to_sieve = m_factorize_upto;
  if (max_prime.has_value())
    max_prime_to_sieve = std::min(max_prime_to_sieve, max_prime.value());

  for (prime_t p = 2; p <= max_prime_to_sieve; ++p) {
    if (m_single_factor.at(p) == p) {
      m_primes.push_back(p);
    }
  }
}

namespace {
inline void factorize_helper(const std::vector<prime_t>& single_factor,
                             Factorization& res, prime_t& v, prime_t v_part) {
  v /= v_part;
  while (v_part != 1) {
    prime_t p = single_factor[v_part];
    v_part /= p;
    // Note that we do not really need the powers of of primes, but we do need
    // to divide by it. We save it for completeness.
    size_t c = 1;
    while (v % p == 0) v /= p, ++c;
    res.push_back({p, c});
  }
}
}  // namespace

void FactorizedSegment::factorize(Factorization& res, prime_t v) {
  const auto& cur = get_value(v);
  res.clear();
  const auto& single_factor = m_base->m_single_factor;
  auto max_prime = m_base->m_max_prime;
  factorize_helper(single_factor, res, v, cur.unique_factors_prod[0]);
  factorize_helper(single_factor, res, v, cur.unique_factors_prod[1]);
  if (v != 1) {
    if (max_prime.has_value() and max_prime.value() < v) return;
    if (v < prime_t(single_factor.size()))
      ASSERT_FATAL(single_factor.at(v) == v);  // Is prime.
    res.emplace_back(v, 1);
  }
}

size_t FactorizedSegment::num_hits_in_sieve(prime_t v) {
  return get_value(v).num_hits;
}

namespace {
using SieveCell = FactorizedSegment::SieveCell;
inline auto& _min(auto& a, auto& b) { return (a < b) ? a : b; }
inline void hit_cell(SieveCell& cell, prime_t threshold, prime_t p) {
  auto& uf = cell.unique_factors_prod;
  auto& min = _min(uf[0], uf[1]);
  if (min <= threshold) min *= p;
  ++cell.num_hits;
}

inline size_t get_start_index(prime_t start, prime_t p) {
  return p - 1 - ((start - 1) % p);
}
}  // namespace

FactorizedSegment::FactorizedSegment(prime_t start,
                                     std::vector<SieveCell> sieve,
                                     std::shared_ptr<FactorizeBase> base)
    : m_start(start), m_sieve(sieve), m_base(base) {}

FactorizedSegment FactorizedSegment::sieve_segment(
    std::shared_ptr<FactorizeBase> base, prime_t start, prime_t end) {
  size_t sieve_size = end - start;
  std::vector<SieveCell> sieve(
      sieve_size, SieveCell{.unique_factors_prod = {1, 1}, .num_hits = 0});
  const auto& primes = base->m_primes;
  const prime_t divisors_bound = base->m_factorize_upto;
  constexpr prime_t small_primes_bound = 1ull << 16;
  // Small primes in this context are primes that we want to sieve in smaller
  // segments (To not re-read the whole array each time).
  {
    // SmallPrimesSinglePass
    auto it =
        std::upper_bound(primes.begin(), primes.end(), small_primes_bound);
    size_t num_small_primes = std::distance(primes.begin(), it);

    std::vector<std::pair<size_t, prime_t>> idx_threshold(num_small_primes);
    for (size_t i = 0; i < num_small_primes; ++i) {
      prime_t p = primes[i];
      idx_threshold[i] = {get_start_index(start, p), divisors_bound / p};
    }
    for (size_t idx_end = 0; idx_end < sieve_size;) {
      idx_end = std::min<size_t>(idx_end + small_primes_bound, sieve_size);
      for (size_t p_index = 0; p_index < num_small_primes; ++p_index) {
        prime_t p = primes[p_index];
        auto& [idx, threshold] = idx_threshold[p_index];
        if (idx >= idx_end) continue;
        // Notice that we update idx_threshold for the next iteration:
        for (; idx < idx_end; idx += p) {
          hit_cell(sieve[idx], threshold, p);
        }
      }
    }
  }
  {
    // BigPrimes
    for (auto p : primes) {
      if (p <= small_primes_bound) continue;
      auto start_index = get_start_index(start, p);
      prime_t threshold = divisors_bound / p;
      for (size_t index = start_index; index < sieve_size; index += p)
        hit_cell(sieve[index], threshold, p);
    }
  }
  return FactorizedSegment(start, sieve, base);
}

const SieveCell& FactorizedSegment::get_value(prime_t v) const {
  ASSERT_FATAL(m_start <= v);
  ASSERT_FATAL(v < m_start + prime_t(m_sieve.size()));
  return m_sieve[v - m_start];
}