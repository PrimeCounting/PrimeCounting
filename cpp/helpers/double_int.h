#pragma once
#include <concepts>
#include <cstdint>
#include <type_traits>

// clang-format off
template <typename T> struct DoubleInt {
  static_assert(sizeof(T) == 0, "T is not one of the possible types...");
};
template <> struct DoubleInt<int32_t>    { using type = int64_t; };
template <> struct DoubleInt<int64_t>    { using type = __int128_t; };
template <> struct DoubleInt<uint32_t>   { using type = uint64_t;};
template <> struct DoubleInt<uint64_t>   { using type = __uint128_t;};

template <typename T> struct HalfInt {
  static_assert(sizeof(T) == 0, "T is not one of the possible types...");
};
template <> struct HalfInt<int64_t>      { using type = int32_t; };
template <> struct HalfInt<__int128_t>   { using type = int64_t; };
template <> struct HalfInt<uint64_t>     { using type = uint32_t;};
template <> struct HalfInt<__uint128_t>  { using type = uint64_t;};
// clang-format on

template <typename T>
using double_int_t = typename DoubleInt<T>::type;

template <typename T>
using half_int_t = typename HalfInt<T>::type;

template <std::integral T>
inline constexpr double_int_t<T> wide_mul(T t1, T t2) {
  return double_int_t<T>(t1) * t2;
}