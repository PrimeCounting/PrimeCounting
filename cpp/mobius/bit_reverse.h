#pragma once

#include <array>
#include <sstream>
#include <vector>

#include "../helpers/assertion.h"

namespace details::bit_reverse {
inline constexpr size_t bit_reverse_naive(size_t i, size_t n) {
  size_t i_rev = 0;
  for (size_t bit = n >> 1; i != 0; i >>= 1, bit >>= 1) {
    if (i & 1) i_rev |= bit;
  }
  return i_rev;
}
inline consteval std::array<uint8_t, 256> reverse_chars() {
  std::array<uint8_t, 256> res;
  for (size_t i = 0; i < 256; ++i) res[i] = bit_reverse_naive(i, 256);
  return res;
}
}  // namespace details::bit_reverse

inline size_t bit_reverse(size_t i, size_t n) {
  ASSERT_FATAL((n & (n - 1)) == 0);
  using namespace details::bit_reverse;
  constexpr auto rev = reverse_chars();
  size_t i_rev = 0;
  for (size_t lg2_n = __builtin_ctzll(n); i != 0; i >>= 8, lg2_n -= 8)
    i_rev |= size_t(rev.at(i & 0xFF)) << lg2_n;
  auto ans = i_rev >> 8;
  return ans;
}

template <typename T>
void bit_reverse_naive(std::vector<T>& v) {
  using namespace details::bit_reverse;
  for (size_t i = 1; i < v.size(); ++i) {
    size_t i_rev = bit_reverse_naive(i, v.size());
    if (i < i_rev) std::swap(v[i], v[i_rev]);
  }
}

template <typename T>
void bit_reverse_impl_small(std::vector<T>& v) {
  const size_t n = v.size();
  ASSERT_FATAL((n & (n - 1)) == 0);
  const auto lg2_n = __builtin_ctzll(n);
  for (size_t i = n - 1, i_rev = n - 1; i > 1; --i) {
    if (i > i_rev) std::swap(v[i], v[i_rev]);
    i_rev ^= (n >> 1);  // We always start with an odd number.
    i--;
    if (i > i_rev) std::swap(v[i], v[i_rev]);
    auto ctz = __builtin_ctzll(i);
    i_rev ^= ((1ull << (ctz + 1)) - 1) << (lg2_n - ctz - 1);
  }
}

namespace details::bit_reverse {
template <size_t B>
consteval auto get_reverse_table() {
  constexpr size_t n = 1ull << B;
  std::array<size_t, n> res;
  for (size_t i = 0; i < n; ++i) res[i] = bit_reverse_naive(i, n);
  return res;
}
}  // namespace details::bit_reverse

template <typename T>
void bit_reverse_impl(std::vector<T>& v) {
  const size_t n = v.size();
  ASSERT_FATAL((n & (n - 1)) == 0);
  const size_t lg2_n = __builtin_ctzll(n);

  constexpr size_t lg2_chunk = 3;
  // For my computer seems that 3 is optimal.
  // TODO: Need to calibrate for different systems.
  if (lg2_n < 2 * lg2_chunk) return bit_reverse_impl_small(v);

  using namespace details::bit_reverse;
  constexpr auto rev = get_reverse_table<lg2_chunk>();
  constexpr size_t chunk = 1 << lg2_chunk;

  size_t chunk_end = n >> lg2_chunk;
  size_t lg2_end = __builtin_ctzll(chunk_end);

  size_t next_mid_rev = 0;
  for (size_t mid = 0, mid_rev = 0; mid < chunk_end;
       mid += chunk, mid_rev = next_mid_rev) {
    size_t bit = chunk_end >> 1;
    while (bit & next_mid_rev) {
      next_mid_rev ^= bit;
      bit >>= 1;
    }
    next_mid_rev ^= bit;

    if (mid > mid_rev) continue;
    if (mid == mid_rev) {
      for (size_t i = 0; i < chunk; ++i) {
        size_t i_rev = rev[i];
        for (size_t j = 0; j < chunk; ++j) {
          size_t j_rev = rev[j];
          size_t idx1 = (i << lg2_end) + mid + j;
          size_t idx2 = (j_rev << lg2_end) + mid + i_rev;
          if (j_rev < i) std::swap(v[idx1], v[idx2]);
        }
      }
    } else {
      for (size_t i = 0; i < chunk; ++i) {
        size_t i_rev = rev[i];
        for (size_t j = 0; j < chunk; ++j) {
          size_t j_rev = rev[j];
          size_t idx1 = (i << lg2_end) + mid + j;
          size_t idx2 = (j_rev << lg2_end) + mid_rev + i_rev;
          std::swap(v[idx1], v[idx2]);
        }
      }
    }
  }
}

template <typename T>
void inplace_bit_reverse(std::vector<T>& v) {
  bit_reverse_impl(v);
}