#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

namespace tqdm {
using clock = std::chrono::steady_clock;
using time_point = clock::time_point;
struct ProgressBar {
  std::string desc;
  std::optional<time_point> begin;
  void set_progress(int prog) {
    if (!begin.has_value()) {
      begin = clock::now();
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                       clock::now() - begin.value())
                       .count();
    std::cerr << desc << ": " << prog << "% [" << elapsed << " (s)] \r";
    if (prog == 100) {
      std::cerr << std::endl;
    }
  }
};

template <typename T>
struct TRange {
  T start, finish, step, idx;
  std::shared_ptr<ProgressBar> pb;
  int64_t last_progress;
  TRange(const TRange&) = default;
  TRange(TRange&&) = default;
  TRange& operator=(const TRange&) = default;
  TRange& operator=(TRange&&) = default;

  TRange(T start, T finish, T step, std::shared_ptr<ProgressBar> pb)
      : start(start),
        finish(finish),
        step(step),
        idx(start),
        pb(pb),
        last_progress(-1) {
    pb->set_progress(0);
  }

  TRange begin() { return *this; }
  TRange end() {
    TRange res(*this);
    res.idx = finish;
    return res;
  }
  bool operator<(const TRange& other) { return *this < other.idx; }
  bool operator!=(const TRange& other) { return *this < other; }

  TRange& operator++() {
    idx += step;
    int32_t new_progress;

    if (*this < finish)
      new_progress =
          std::floor(abs_diff(idx, start) * 100. / abs_diff(finish, start));
    else
      new_progress = 100;

    if (last_progress != new_progress) {
      pb->set_progress(new_progress);
      last_progress = new_progress;
    }
    return *this;
  }
  T operator*() const { return idx; }

 private:
  static T abs_diff(T a, T b) { return a < b ? b - a : a - b; }
  bool operator<(T other_idx) {
    if (step > 0)
      return idx < other_idx;
    else
      return idx > other_idx;
  }
};

inline auto get_pb() { return std::make_shared<ProgressBar>(); }

template <typename T>
inline auto range(T n) {
  return TRange<T>(0, n, 1, get_pb());
}
template <typename T>
inline auto range(T start, T finish) {
  return TRange<T>(start, finish, 1, get_pb());
}
template <typename T>
inline auto range(T start, T finish, T step) {
  return TRange<T>(start, finish, step, get_pb());
}

template <typename T, typename... Args>
auto title_range(std::string desc, Args... start_end_step) {
  auto res = range<T>(start_end_step...);
  res.pb->desc = desc;
  res.pb->set_progress(0);
  return res;
}

struct Title {
  Title(std::string desc) : tq(title_range<size_t>(desc, 1)) {}
  ~Title() { ++tq; }
  TRange<size_t> tq;
};

}  // namespace tqdm