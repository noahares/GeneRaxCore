#pragma once

#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>

/**
 *  Return the indices of the input vector sorted in descending order
 */
template <typename T>
std::vector<size_t> sort_indices_descending(const std::vector<T> &v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
  return idx;
}

struct Counter {
  struct value_type {
    template <typename T> value_type(const T &) {}
  };
  void push_back(const value_type &) { ++count; }
  size_t count = 0;
};

/**
 *  Return the number of elements present in both sorted containers
 */
template <typename T1, typename T2>
size_t intersection_size(const T1 &s1, const T2 &s2) {
  Counter c;
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(c));
  return c.count;
}
