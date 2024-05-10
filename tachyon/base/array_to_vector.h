#ifndef TACHYON_BASE_ARRAY_TO_VECTOR_H_
#define TACHYON_BASE_ARRAY_TO_VECTOR_H_

#include <stddef.h>

#include <vector>

namespace tachyon::base {

template <typename T>
const std::vector<T> ArrayToVector(const T (&arr)[0]) {
  return {};
}

template <typename T, size_t N>
const std::vector<T> ArrayToVector(const T (&arr)[N]) {
  return std::vector<T>(std::begin(arr), std::end(arr));
}

template <typename T>
std::vector<std::vector<T>> Array2DToVector2D(const T (&arr)[0][0]) {
  return {};
}

template <typename T, size_t A, size_t B>
std::vector<std::vector<T>> Array2DToVector2D(const T (&arr)[A][B]) {
  std::vector<std::vector<T>> vec;
  vec.reserve(A);
  for (const auto& inner_array : arr) {
    vec.emplace_back(std::begin(inner_array), std::end(inner_array));
  }
  return vec;
}

}  // namespace tachyon::base

#endif  // TACHYON_BASE_ARRAY_TO_VECTOR_H_
