#include <iostream>
#include <sstream>

#include "count_primes/count_primes.h"

int main(int argc, char* argv[]) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: countprimes UPTO [MEMORY_TRADEOFF]" << std::endl;
    return -1;
  }
  long double upto_;
  std::stringstream s1(argv[1]);
  s1 >> upto_;
  prime_t upto = std::round(upto_);
  if (s1.fail() || upto < 2 || upto_ != upto) {
    std::cerr << "Could not parse first argument." << std::endl;
    return -1;
  }
  std::cout << upto << std::endl;
  double memory_tradeoff = 1.;
  if (argc == 3) {
    std::stringstream s2(argv[2]);
    s2 >> memory_tradeoff;
    if (s2.fail() || memory_tradeoff < 0) {
      std::cout << "Could not parse second argument." << std::endl;
      return -1;
    }
  }
  double lg2_prec = 1. / std::sqrt(upto) * memory_tradeoff;
  prime_t max_prime_to_use = std::ceil(std::sqrt(upto));

  auto computed_num_primes = count_primes(upto, lg2_prec, max_prime_to_use);
  std::cout << "Num primes up to " << upto << ":" << std::endl
            << "\t" << computed_num_primes << std::endl;
  return 0;
}