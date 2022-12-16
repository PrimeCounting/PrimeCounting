#pragma once
#include <concepts>
#include <iostream>

#include "../constants.h"
#include "double_int.h"
#include "math.h"
#include "types.h"

namespace details {
template <std::unsigned_integral T, T MOD>
static consteval T compute_r() {
  T res = 1;
  while (res * MOD != 1) {
    // res*MOD = 1+2a
    // res*MOD*(1-2a) = 1-4a^2
    res *= 2 - res * MOD;
  }
  return res;
}
template <std::unsigned_integral T, T MOD>
struct ModInt_impl {
  using SignedT = std::make_signed_t<T>;
  static_assert(SignedT(MOD) * 2 == wide_mul<T>(MOD, 2),
                "Need to be able to contain 2 * MOD");
  static constexpr T mod = MOD;
  static_assert(is_prime(MOD));
  using Base = ModInt_impl<T, MOD>;

  constexpr ModInt_impl() : value(0){};

  constexpr ModInt_impl(const Base &) = default;
  constexpr ModInt_impl(Base &&) = default;

  template <std::integral U>
  constexpr ModInt_impl(U v)
      : value(montgomery_mod(wide_mul(my_mod(v), two_power_2num_bits))) {}

  constexpr Base &operator=(const Base &) = default;
  constexpr Base &operator=(Base &&) = default;

  constexpr bool operator==(const Base &other) const {
    return normalized() == other.normalized();
  }
  constexpr Base &operator*=(const Base &other) {
    value = montgomery_mod(wide_mul<T>(value, other.value));
    return *this;
  }
  constexpr Base &operator+=(const Base &other) {
    auto new_v = SignedT(value) - 2 * SignedT(MOD) + SignedT(other.value);
    // We use signed here because comparison with 0 is faster.
    // We use signed here because comparison with 0 is faster.
    if (new_v < 0) new_v += 2 * MOD;
    value = new_v;
    return *this;
  }
  constexpr Base &operator-=(const Base &other) {
    auto new_v = SignedT(value) - SignedT(other.value);
    if (new_v < 0) new_v += 2 * MOD;
    value = new_v;
    return *this;
  }
  constexpr Base &operator/=(const Base &o) { return *this *= o.inverse(); }

  constexpr Base operator*(const Base &o) const { return Base(*this) *= o; }
  constexpr Base operator+(const Base &o) const { return Base(*this) += o; }
  constexpr Base operator-(const Base &o) const { return Base(*this) -= o; }
  constexpr Base operator/(const Base &o) const { return Base(*this) /= o; }

  constexpr Base &operator++() { return *this += 1; }
  constexpr Base &operator--() { return *this -= 1; }
  constexpr Base operator-() const { return Base(*this).neg(); }
  constexpr Base &neg() {
    value = 2 * MOD - value;
    return *this;
  }
  constexpr Base pow(size_t power) const {
    Base res = 1;
    for (Base base = *this; power != 0; power /= 2, base *= base) {
      if (power & 1) res *= base;
    }
    return res;
  }

  constexpr Base inverse() const { return this->pow(MOD - 2); }

  // Return the value modulo MOD ([0, MOD))
  constexpr T get() const {
    T ans = montgomery_mod(value);
    // montgomery_mod might return value in [0,2*MOD):
    return ans >= MOD ? ans - MOD : ans;
  }

  friend std::ostream &operator<<(std::ostream &out, const Base &v) {
    return out << "ModInt<" << MOD << ">{" << v.get() << "}";
  }

  static constexpr T get_mod() { return MOD; }

 private:
  static constexpr size_t num_bits = sizeof(T) * 8;
  static constexpr T two_power_2num_bits = pow_mod(T(2), num_bits * 2, MOD);
  static constexpr T r = compute_r<T, MOD>();
  static constexpr inline T montgomery_mod(double_int_t<T> v) {
    // Assume v = a * 2 ^ 2num_bits
    // Returns a * 2 ^ num_bits % MOD
    T ans = (v >> num_bits) +
            (T(v) != 0);  // Whether there is a carry frow lower bits.
    ans += (wide_mul<T>(T(-r) * T(v), MOD)) >> num_bits;
    return ans;
  }
  constexpr T normalized() const { return value >= MOD ? value - MOD : value; }

  // Return value in range [0,2MOD)
  template <std::signed_integral U>
  constexpr T my_mod(U v) {
    return U(v % SignedT(MOD)) + SignedT(MOD);
  }
  template <std::unsigned_integral U>
  constexpr T my_mod(U v) {
    return v % MOD;
  }
  // equivalent to real_value * 2^num_bits
  // in range [0,2*MOD)
  T value;
};
}  // namespace details
template <__uint128_t MOD>
using ModInt128 = details::ModInt_impl<__uint128_t, MOD>;
template <uint64_t MOD>
using ModInt64 = details::ModInt_impl<uint64_t, MOD>;
template <uint32_t MOD>
using ModInt32 = details::ModInt_impl<uint32_t, MOD>;

template <std::unsigned_integral T, T MOD>
struct std::is_integral<details::ModInt_impl<T, MOD>> : std::true_type {};

using mint = ModInt64<MOD>;
