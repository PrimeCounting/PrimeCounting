#include "count_primes.h"

#include "../factorize_range/factorize_range.h"
#include "../helpers/cell.h"
#include "../helpers/double_int.h"
#include "../helpers/indicators.h"
#include "../helpers/sieve_primes.h"
#include "../mobius/mobius_using_newton.h"
#include "logarithmic_integral.h"

mint count_primes_with_errors(prime_t upto, double lg2_prec,
                              prime_t max_prime_to_use) {
  size_t num_small_primes = get_primes_by_sieve(max_prime_to_use).size();
  // Mobius only of numbers with factors are up to max_prime_to_use.
  auto mobius = get_mobius_using_newton(upto, lg2_prec,
                                        /*max_prime=*/max_prime_to_use);

  auto get_cumsum_all_numbers = [lg2_prec](size_t cell) {
    return get_cell_end(cell, lg2_prec);
  };
  /* We want to compute:
   * large_primes = conv(all_numbers, mobius);
   * ans = sum(large_primes[:max_cell+1]);
   * This is equivalent to multiplying each index `i` in mobius by the cumsum of
   * all_numbers upto index `max_cell - i`.
   */
  size_t max_cell = get_cell(upto, lg2_prec);
  mint estimated_num_large_primes = 0;  // Larger than `max_prime_to_use`
  {
    // Bone*Mobius
    size_t i = 0;
    size_t last_cumsum = 2 * upto;
    {
      // Bone*Mobius_sparse_part
      constexpr size_t dense_bone_cell_size_bound = 10;
      for (; i <= max_cell; ++i) {
        if (mobius[i] != 0) {
          // mobius[:mid_cell] is mostly zeros.
          size_t cur_cumsum = get_cumsum_all_numbers(max_cell - i);
          estimated_num_large_primes += mobius[i] * cur_cumsum;
          size_t cur_cell_size = last_cumsum - cur_cumsum;
          last_cumsum = cur_cumsum;
          if (cur_cell_size < dense_bone_cell_size_bound) {
            // We are in the dense part, and \bone is in the sparse part!
            ++i;
            break;
          }
        }
      }
    }
    {
      // Bone*Mobius_dense_part
      size_t next_nnz = last_cumsum;
      size_t next_cell = get_cell(next_nnz, lg2_prec);
      for (; i <= max_cell; ++i) {
        size_t cur_cell = max_cell - i;
        while (cur_cell < next_cell) {
          --next_nnz;
          next_cell = get_cell(next_nnz, lg2_prec);
        }
        estimated_num_large_primes += mobius[i] * next_nnz;
      }
    }
  }
  // `estimated_num_large_primes` also includes 1.
  return estimated_num_large_primes + num_small_primes - 1;
}

namespace {
prime_t get_closest(prime_t v, std::array<prime_t, 3> candidates) {
  prime_t ans = candidates[0];
  for (auto i : candidates) {
    if (std::abs(i - v) < std::abs(ans - v)) ans = i;
  }
  return ans;
}

prime_t lift_to_integer_using_li(mint prime_count, prime_t upto) {
  // Assuming the Riemann Hypothesis, we can deduce the top bits of the result
  // using the logarithm integral.
  prime_t estimation = std::round(li(upto));
  prime_t v = prime_count.get();
  prime_t mod = mint::get_mod();
  prime_t est_mod = estimation - (estimation % mod);
  return get_closest(estimation,
                     {est_mod - mod + v, est_mod + v, est_mod + mod + v});
}
}  // namespace

prime_t count_primes(prime_t upto, double lg2_prec, prime_t max_prime_to_use) {
  mint ans = count_primes_with_errors(upto, lg2_prec, max_prime_to_use);
  ans -= error_correction(upto, lg2_prec, max_prime_to_use);
  return lift_to_integer_using_li(ans, upto);
}