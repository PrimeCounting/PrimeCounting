
#include "count_primes.h"

#include <gtest/gtest.h>

#include "../helpers/sieve_primes.h"

void test_size(prime_t upto) {
  prime_t num_primes = get_primes_by_sieve(upto).size();
  EXPECT_EQ(num_primes, count_primes(upto));
}

TEST(correctness_test, high_precision) {
  /**
   * This test is intendant to verify that given high enough precision the
   * answer should be correct even without error corrections.
   */
  constexpr prime_t upto = 900;
  constexpr prime_t use_primes_upto = 30;
  constexpr double lg2_prec = 1. / upto;

  size_t num_primes = get_primes_by_sieve(upto).size();
  auto expected = count_primes_with_errors(upto, lg2_prec, use_primes_upto);
  EXPECT_EQ(num_primes, expected);
}

TEST(correctness_test, test_upto_2pow25) {
  for (prime_t lg2_upto = 1; lg2_upto <= 25; ++lg2_upto)
    test_size(1ull << lg2_upto);
}

TEST(correctness_test, test_upto_1000) {
  for (prime_t upto = 2; upto <= 1000; ++upto) test_size(upto);
}

TEST(correctness_test, error_correction_small) { test_size(10); }

__int128_t parse_string(std::string v) {
  __int128_t value = 0;
  for (char c : v) {
    ASSERT_FATAL('0' <= c and c <= '9');
    value *= 10;
    value += c - '0';
  }
  return value;
}
std::map<size_t, __int128_t> ans_from_web(
    {{10, parse_string("455052511")},
     {11, parse_string("4118054813")},
     {12, parse_string("37607912018")},
     {13, parse_string("346065536839")},
     {14, parse_string("3204941750802")},
     {15, parse_string("29844570422669")},
     {16, parse_string("279238341033925")},
     {17, parse_string("2623557157654233")},
     {18, parse_string("24739954287740860")},
     {19, parse_string("234057667276344607")},
     {20, parse_string("2220819602560918840")},
     {21, parse_string("21127269486018731928")},
     {22, parse_string("201467286689315906290")}});

TEST(correctness_test, count_primes) {
  constexpr size_t power = 14;
  prime_t upto = pow10(power);

  double lg2_prec = 1. / std::sqrt(upto) * 5;
  prime_t use_primes_upto = std::floor(std::sqrt(upto));

  auto computed_num_primes = count_primes(upto, lg2_prec, use_primes_upto);
  if (power < 10) {
    EXPECT_EQ(get_primes_by_sieve(upto).size(), computed_num_primes);
  } else {
    if (mint(ans_from_web[power]) != computed_num_primes) {
      EXPECT_EQ(ans_from_web[power], computed_num_primes);
    }
  }
}