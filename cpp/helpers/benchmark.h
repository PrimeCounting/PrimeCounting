#pragma once
#include <algorithm>
#include <chrono>

inline void benchmark_code(auto&& func, double max_time = 20) {
  double total_time = 0;
  double min_time = 1e300;
  for (size_t i = 0; i < 30 && total_time < max_time; ++i) {
    auto start_time = std::chrono::system_clock::now();
    auto res = func();
    if (res == 123342421231432ull)
      std::cerr << res << std::endl;  // should disable optimizations :)
    auto end_time = std::chrono::system_clock::now();
    auto t = std::chrono::duration<double>(end_time - start_time).count();
    total_time += t;
    min_time = std::min(min_time, t);
  }
  std::cout << "Min run of code took: " << min_time << std::endl;
}
