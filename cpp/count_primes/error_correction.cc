#include "error_correction.h"

#include <map>

#include "../factorize_range/factorize_range.h"
#include "../helpers/cell.h"
#include "../helpers/indicators.h"
#include "../helpers/math.h"

namespace {

struct Factor {
  prime_t p;
  size_t cell;
};

struct CorrectionRecursion {
  CorrectionRecursion(prime_t v, const std::vector<Factor>& factor_cell_pairs,
                      double lg2_prec, size_t max_cell)
      : m_v(v),
        m_factor_cell_pairs(factor_cell_pairs),
        m_lg2_prec(lg2_prec),
        m_max_cell(max_cell) {}

  prime_t compute() const {
    size_t cell = 0, mob = 1;
    prime_t prime_mult = 1;
    for (auto [p, c] : m_factor_cell_pairs) {
      cell += c;
      prime_mult *= p;
      mob *= -1;
    }
    return recursion(prime_mult, cell, mob, 0);
  }

 private:
  prime_t recursion(prime_t curr, size_t curr_cell, int32_t curr_mobius,
                    size_t factor_index) const {
    size_t rem_cell = get_cell(m_v / curr, m_lg2_prec);
    size_t affected_cell = rem_cell + curr_cell;
    if (affected_cell > m_max_cell) return 0;
    prime_t divisors_error = curr_mobius;
    for (size_t i = factor_index; i < m_factor_cell_pairs.size(); ++i) {
      auto [prime, cell] = m_factor_cell_pairs[i];
      auto new_curr = curr / prime;
      auto new_cell = curr_cell - cell;
      divisors_error += recursion(new_curr, new_cell, curr_mobius * -1, i + 1);
    }
    return divisors_error;
  }

  prime_t m_v;
  const std::vector<Factor>& m_factor_cell_pairs;
  double m_lg2_prec;
  size_t m_max_cell;
};

prime_t error_correction(prime_t v, size_t max_cell, double lg2_prec,
                         const Factorization& factors) {
  static thread_local std::vector<Factor> factor_cell_vec;
  factor_cell_vec.clear();
  for (auto& p_c : factors) {
    auto p = p_c.first;
    factor_cell_vec.push_back({p, get_cell(p, lg2_prec)});
  }
  return CorrectionRecursion(v, factor_cell_vec, lg2_prec, max_cell).compute();
}

struct CachedGetCell {
  CachedGetCell(prime_t start, prime_t end, double lg2_prec) : m_end(end) {
    ASSERT_FATAL(start < end);
    prime_t cur = start;
    prime_t last_cur;
    do {
      last_cur = cur;
      size_t cur_cell = get_cell(cur, lg2_prec);
      thresholds[cur] = cur_cell;
      cur = get_cell_end(cur_cell, lg2_prec) + 1;
      ASSERT_FATAL(thresholds.size() < 100u);
    } while (last_cur < end);
  }

  size_t operator()(prime_t v) const {
    ASSERT_FATAL(v < m_end);
    auto it = thresholds.lower_bound(v);  // Notice reversed order.
    ASSERT_FATAL(it != thresholds.end());
    return it->second;
  };
  std::map<prime_t, size_t, std::greater<prime_t>> thresholds;
  prime_t m_end;
};
}  // namespace

prime_t error_correction(prime_t upto, double lg2_prec,
                         prime_t max_prime_to_use) {
  size_t max_cell = get_cell(upto, lg2_prec);
  prime_t max_value_to_check = get_max_value_to_check(max_cell, lg2_prec);
  CachedGetCell cached_get_cell(upto + 1, max_value_to_check + 1, lg2_prec);

  auto should_factor = [&](prime_t v, size_t num_hits) -> bool {
    return max_cell + num_hits >= cached_get_cell(v);
  };
  struct Error {
    prime_t amount;
    size_t num_errors;
    Error& operator+=(const Error& other) {
      amount += other.amount;
      num_errors += other.num_errors;
      return *this;
    }
  };
  auto handle_error = [&](prime_t v, Factorization& factors) -> Error {
    return Error{.amount = error_correction(v, max_cell, lg2_prec, factors),
                 .num_errors = 1u};
  };
  auto res = handle_range(upto + 1, max_value_to_check + 1, should_factor,
                          handle_error, "Error correction", max_prime_to_use);
  return res.amount;
}  // namespace
