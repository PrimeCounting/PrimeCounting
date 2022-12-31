#include "logarithmic_integral.h"

#include <cmath>

#include "../helpers/types.h"

double li(prime_t upto) {
  if (upto <= 2) return 1;
  if (upto <= 3) return 2;
  const long double x = upto;
  const long double lnx = std::log(x);
  const long double lnlnx = std::log(lnx);

  long double ans = lnlnx;
  // Euler-Mascheroni constant
  ans += 0.577215664901532l;
  // sqrt(x) * (-1)^(n-1) * (lnx)^n / (n!2^(n-1))
  long double err_term_1 = -2 * std::sqrt(x);
  // sum(k=1...(n-1)/2) 1/(2k+1)
  long double err_term_2 = 0;
  for (size_t n = 1;; ++n) {
    err_term_1 *= -lnx / (n * 2);
    if (n & 1) {
      err_term_2 += 1.0L / n;
    }
    auto err = err_term_1 * err_term_2;
    ans += err;
    if (std::abs(err) < 1e-7) {
      break;
    }
  }
  return ans;
}
