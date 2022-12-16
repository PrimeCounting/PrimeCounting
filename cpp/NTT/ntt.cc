#include "ntt.h"

#include <bits/stdc++.h>

#include "../helpers/assertion.h"
#include "../helpers/indicators.h"
#include "../mobius/bit_reverse.h"

namespace {

template <typename mint>
consteval std::pair<mint, size_t> get_root() {
  auto mod = mint::get_mod();
  auto ord = mod - 1;
  decltype(ord) part2 = 1;
  size_t ord2 = 0;
  while ((ord % 2) == 0) {
    part2 *= 2;
    ord /= 2;
    ++ord2;
  }

  for (mint res = 2;; ++res) {
    auto res1 = res.pow(ord);
    if (res1.pow(part2 / 2) != 1) {
      return {res1, ord2};
    }
  }
}

inline void butterfly(auto& vec, size_t block_start, size_t block_size,
                      size_t chunk_size, const auto& w_powers) {
  for (size_t i = 0; i < chunk_size; i++) {
    auto u = vec[block_start + i];
    auto v = vec[block_start + i + block_size / 2] * w_powers[i];
    vec[block_start + i] = u + v;
    vec[block_start + i + block_size / 2] = u - v;
  }
}

template <typename T, size_t N>
inline void generate_w_powers(std::array<T, N>& w_powers, T w,
                              size_t max_power) {
  constexpr size_t batch = 1ull << 3;
  static_assert(batch <= N);

  w_powers[0] = 1;
  for (size_t i = 1; i < batch; ++i) w_powers[i] = w * w_powers[i - 1];
  const auto jump = w.pow(batch);
  for (size_t i = batch; i < max_power; i += batch)
    for (size_t j = 0; j < batch; ++j)
      w_powers[i + j] = jump * w_powers[i - batch + j];
}
}  // namespace

void details::ntt_impl(std::vector<mint>& vec, bool inverse,
                       std::optional<std::string> desc) {
  constexpr auto root_info = get_root<mint>();

  auto [root, ord2] = root_info;
  while (vec.size() != (1ull << ord2)) {
    root = root.pow(2);
    --ord2;
    ASSERT_FATAL(ord2 != 0);  // Not a power of two!
  }
  if (inverse) root = root.inverse();

  inplace_bit_reverse(vec);
  size_t fft_size = vec.size();

  std::optional<tqdm::TRange<size_t>> tq;
  if (desc.has_value()) {
    tq = tqdm::title_range<size_t>(desc.value(), ord2 - 1);
  }
  // We compute the w_powers (the coefficient for the butterfly) in chunks.
  constexpr size_t in_block_chunk_size = 1ull << 3;
  // We reuse the same chunk across multiple blocks to reduce computation.
  constexpr size_t num_parallel_blocks = 1ull << 2;

  std::array<mint, in_block_chunk_size> w_powers;
  for (size_t block_size = 2; block_size <= fft_size; block_size *= 2) {
    const auto cur_w = root.pow(fft_size / block_size);
    if (block_size / 2 <= in_block_chunk_size) {
      // SmallBlock
      generate_w_powers(w_powers, cur_w, block_size / 2);
      for (size_t block_start = 0; block_start < fft_size;
           block_start += block_size)
        butterfly(vec, block_start, block_size, block_size / 2, w_powers);
    } else {
      // LargeBlock
      const mint jump = cur_w.pow(in_block_chunk_size);
      const size_t parallel_blocks_size =
          std::min(fft_size, block_size * num_parallel_blocks);
      for (size_t idx = 0; idx < fft_size; idx += parallel_blocks_size) {
        generate_w_powers(w_powers, cur_w, in_block_chunk_size);
        for (size_t in_block_idx = 0; in_block_idx < block_size / 2;
             in_block_idx += in_block_chunk_size) {
          for (size_t j = 0; j < parallel_blocks_size; j += block_size) {
            size_t start = idx + j + in_block_idx;
            butterfly(vec, start, block_size, in_block_chunk_size, w_powers);
          }
          for (size_t i = 0; i < in_block_chunk_size; ++i) w_powers[i] *= jump;
        }
      }
    }
    if (tq.has_value()) {
      ++tq.value();
    }
  }
  if (inverse) {
    mint v = mint(fft_size).inverse();
    for (auto& i : vec) i *= v;
  }
}