#pragma once
#include <stdexcept>

// clang-format off
#define ASSERT_FATAL(cond)   \
  if (!(cond))               \
  throw std::runtime_error(  \
    std::string(#cond) + "\n" + __FILE__ + ":" + std::to_string(__LINE__))
// clang-format on