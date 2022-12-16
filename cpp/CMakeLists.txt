cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(CountPrimes VERSION 1.0)


# CPP flags:
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Werror -march=native -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


add_library(
     count_primes
     all_numbers/all_numbers.cc
     count_primes/count_primes.cc
     count_primes/error_correction.cc
     count_primes/logarithmic_integral.cc
     factorize_range/factorize_range.cc
     helpers/cell.cc
     helpers/math.cc
     helpers/sieve_primes.cc
     mobius/mobius_using_newton.cc
     NTT/ntt.cc
)

add_executable(countprimes "countprimes.cc")
target_link_libraries(countprimes count_primes)

option(COMPILE_TESTS "Whether to compile the tests" OFF)

# TESTS:
if (COMPILE_TESTS)
  find_package(GTest REQUIRED)
  include(GoogleTest)
  file(GLOB_RECURSE tests_files "*_test.cc")
  add_executable(tests ${tests_files})
  target_link_libraries(tests count_primes)
  target_link_libraries(tests GTest::GTest GTest::Main)
  gtest_discover_tests(tests)
endif()