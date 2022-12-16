#pragma once
#include <stdint.h>

#include <memory>
#include <vector>

#include "../helpers/assertion.h"
#include "../helpers/mod_int.h"

namespace details {
void ntt_impl(std::vector<mint>& v, bool inverse,
              std::optional<std::string> desc = std::nullopt);
}  // namespace details

template <typename... Args>
void ntt(std::vector<mint>& v, Args... args) {
  return details::ntt_impl(v, /*inverse=*/false, args...);
}

template <typename... Args>
void intt(std::vector<mint>& v, Args... args) {
  return details::ntt_impl(v, /*inverse=*/true, args...);
}